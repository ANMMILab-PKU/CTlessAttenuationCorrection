function [FI, alpha_rec] = calc_fisher_information_SCS(atten_map, activity_map, file_em, ...
    file_tr_202, file_bs_202, file_tr_307, file_bs_307, params_SCS, psi_geom)
% CALC_FISHER_INFORMATION_SCS
% -------------------------------------------------------------------------
% 计算 Emission 和 Transmission 数据源的近似 Fisher 信息量，
% 用于自动平衡各数据源的权重 alpha。
%
% 理论依据:
%   对于 Poisson 模型 y ~ Poisson(w)，Fisher 信息量:
%     FI_j = sum_i A_ij^2 * u_i^2 / w_i
%   利用可分离上界 A_ij^2 <= A_ij * L_i，近似:
%     FI_j ≈ [A^T(u^2/w * L)]_j
%   总 Fisher 信息: sum_j FI_j ≈ sum_i (u_i^2 / w_i) * L_i^2
%
%   对于 Transmission，需乘以能量映射链式法则导数:
%     FI_j^TX = eta'_j^2 * [A^T(u^2/w * L)]_j
%   全局近似: sum_j FI_j^TX ≈ mean(eta'^2) * sum_i (u_i^2/w_i * L_i^2)
%
% 输出:
%   FI:        结构体，包含各数据源的 Fisher 信息量
%   alpha_rec: 结构体，包含推荐的平衡权重
%
% 用法:
%   在主迭代循环之前调用一次:
%   [FI, alpha_rec] = calc_fisher_information_SCS(atten_map, activity_map, ...);
%   params_SCS.alpha_EM  = alpha_rec.alpha_EM;
%   params_SCS.alpha_202 = alpha_rec.alpha_202 * gamma;  % gamma 为用户调节因子
%   params_SCS.alpha_307 = alpha_rec.alpha_307 * gamma;
% -------------------------------------------------------------------------

% 解包几何参数
libname    = psi_geom.libname;
img_dim    = psi_geom.img_dim;
img_origin = psi_geom.img_origin;
voxel_size = psi_geom.voxel_size;
chunk_size = params_SCS.chunk_size;
Nvox       = prod(img_dim);
all_1_map  = ones(Nvox, 1, 'single');

FI_em    = double(0);
FI_tx202 = double(0);
FI_tx307 = double(0);

fprintf('\n=== Computing Fisher Information for Weight Balancing ===\n');

% =====================================================================
% 1. Dense Emission Pass (所有几何 LOR 的基线 Fisher 信息)
% =====================================================================
if params_SCS.alpha_EM > 0
    fprintf('  [FI] Processing Dense Emission Geometry...\n');
    for i_part = 1 : psi_geom.num_parts
        sino_id = Get_valid_lor_calculate_sens_single_layer( ...
            psi_geom.ring_num, psi_geom.crystal_per_ring, ...
            psi_geom.min_crystal_difference, psi_geom.max_ring_difference, ...
            psi_geom.num_parts, i_part);

        xs = single(psi_geom.crystal_lut(sino_id(:,1), 1:3)'); xs = xs(:);
        xe = single(psi_geom.crystal_lut(sino_id(:,2), 1:3)'); xe = xe(:);
        nlors = size(sino_id, 1);
        s1 = ones(nlors, 1, 'single');

        proj_mu  = proj_forw_parallelproj_cuda(libname, xs, xe, atten_map,    img_origin, voxel_size, s1, nlors, img_dim);
        proj_act = proj_forw_parallelproj_cuda(libname, xs, xe, activity_map,  img_origin, voxel_size, s1, nlors, img_dim);
        L_path   = proj_forw_parallelproj_cuda(libname, xs, xe, all_1_map,     img_origin, voxel_size, s1, nlors, img_dim);

        nf = NF_generate_data(sino_id - 1, psi_geom.crystal_lut, ...
            psi_geom.PET_geom, psi_geom.NF_ax_data, psi_geom.NF_tr_data);

        u_i = exp(-proj_mu) .* proj_act ./ nf;
        u_i(isnan(u_i) | isinf(u_i)) = 0;

        % Dense pass: w ≈ u (no randoms/scatter), FI_i ≈ u_i * L_i^2
        FI_em = FI_em + sum(double(u_i .* L_path.^2));
    end
end

% =====================================================================
% 2. Sparse Emission Correction (y > 0 的 LOR 修正)
% =====================================================================
if params_SCS.alpha_EM > 0
    fprintf('  [FI] Processing Sparse Emission Histogram...\n');
    fid_em = fopen(file_em, 'rb');
    if fid_em == -1, warning('Cannot open emission file, skipping sparse FI.');
    else
        while ~feof(fid_em)
            d_em = fread(fid_em, [3, chunk_size], 'uint32');
            if isempty(d_em), break; end
            n_em = size(d_em, 2);

            id1 = d_em(2,:) + 1;  id2 = d_em(3,:) + 1;
            xs = single(psi_geom.crystal_lut(id1, 1:3)'); xs = xs(:);
            xe = single(psi_geom.crystal_lut(id2, 1:3)'); xe = xe(:);
            s1 = ones(n_em, 1, 'single');

            proj_mu  = proj_forw_parallelproj_cuda(libname, xs, xe, atten_map,   img_origin, voxel_size, s1, n_em, img_dim);
            proj_act = proj_forw_parallelproj_cuda(libname, xs, xe, activity_map, img_origin, voxel_size, s1, n_em, img_dim);
            L_path   = proj_forw_parallelproj_cuda(libname, xs, xe, all_1_map,    img_origin, voxel_size, s1, n_em, img_dim);

            nf = NF_generate_data([d_em(2,:)', d_em(3,:)'], psi_geom.crystal_lut, ...
                psi_geom.PET_geom, psi_geom.NF_ax_data, psi_geom.NF_tr_data);
            [r_em, ~, ~] = get_additive_terms(id1, id2, n_em);

            u_i = exp(-proj_mu) .* proj_act ./ nf;
            u_i(isnan(u_i) | isinf(u_i)) = 0;
            w_i = max(u_i + single(r_em), 1e-10);

            % 修正: 真实 FI 为 u^2/w * L^2，基线为 u * L^2
            % 修正量 = (u^2/w - u) * L^2 = u * L^2 * (u/w - 1) < 0
            FI_em = FI_em + sum(double((u_i.^2 ./ w_i - u_i) .* L_path.^2));
        end
        fclose(fid_em);
    end
end

% =====================================================================
% 3. Transmission 202 keV Fisher Information
% =====================================================================
if params_SCS.tx_mode == 1 || params_SCS.tx_mode == 3
    fprintf('  [FI] Processing Transmission 202 keV...\n');
    [mu_low_202, eta_202] = map_mu511_to_low_energy_3segLM(atten_map, 202);
    mean_eta_sq_202 = mean(double(eta_202(:)).^2);

    FI_tx202 = compute_tx_fisher( ...
        file_tr_202, file_bs_202, params_SCS.file_tr_202_add, ...
        mu_low_202, mean_eta_sq_202, params_SCS, psi_geom);
end

% =====================================================================
% 4. Transmission 307 keV Fisher Information
% =====================================================================
if params_SCS.tx_mode == 2 || params_SCS.tx_mode == 3
    fprintf('  [FI] Processing Transmission 307 keV...\n');
    [mu_low_307, eta_307] = map_mu511_to_low_energy_3segLM(atten_map, 307);
    mean_eta_sq_307 = mean(double(eta_307(:)).^2);

    FI_tx307 = compute_tx_fisher( ...
        file_tr_307, file_bs_307, params_SCS.file_tr_307_add, ...
        mu_low_307, mean_eta_sq_307, params_SCS, psi_geom);
end

% =====================================================================
% 5. 汇总输出与推荐权重
% =====================================================================
FI.emission = FI_em;
FI.tx_202   = FI_tx202;
FI.tx_307   = FI_tx307;
FI.total    = FI_em + FI_tx202 + FI_tx307;

% 推荐权重: 使各数据源的 Fisher 信息贡献相当
% alpha_EM 固定为 1.0，TX 权重按 FI 比值调整
alpha_rec.alpha_EM  = 1.0;
alpha_rec.alpha_202 = FI_em / max(FI_tx202, eps);
alpha_rec.alpha_307 = FI_em / max(FI_tx307, eps);

fprintf('\n=== Fisher Information Results ===\n');
fprintf('  FI_emission  = %.6e\n', FI_em);
fprintf('  FI_TX_202    = %.6e\n', FI_tx202);
fprintf('  FI_TX_307    = %.6e\n', FI_tx307);
fprintf('  FI ratio (EM/TX202) = %.4f\n', FI_em / max(FI_tx202, eps));
fprintf('  FI ratio (EM/TX307) = %.4f\n', FI_em / max(FI_tx307, eps));
fprintf('  Recommended balanced weights:\n');
fprintf('    alpha_EM  = %.4f\n', alpha_rec.alpha_EM);
fprintf('    alpha_202 = %.4f\n', alpha_rec.alpha_202);
fprintf('    alpha_307 = %.4f\n', alpha_rec.alpha_307);
fprintf('  Note: multiply TX weights by gamma (1~3) to favor TX data.\n');
fprintf('=========================================\n\n');
end

% =========================================================================
% 内部辅助函数: 计算单个 TX 数据源的 Fisher 信息量
% =========================================================================
function FI_tx = compute_tx_fisher(f_tr, f_bs, f_add, mu_low_img, mean_eta_sq, params, geom)
FI_tx = double(0);
Nvox = prod(geom.img_dim);
all_1_map = ones(Nvox, 1, 'single');

fid_tr  = fopen(f_tr,  'rb');
fid_bs  = fopen(f_bs,  'rb');
fid_add = fopen(f_add, 'rb');

if fid_tr == -1 || fid_bs == -1 || fid_add == -1
    warning('Cannot open TX files. Returning FI = 0.');
    if fid_tr  ~= -1, fclose(fid_tr);  end
    if fid_bs  ~= -1, fclose(fid_bs);  end
    if fid_add ~= -1, fclose(fid_add); end
    return;
end

while ~feof(fid_tr)
    d_tr  = fread(fid_tr,  [3, params.chunk_size], 'uint32');
    d_bs  = fread(fid_bs,  [3, params.chunk_size], 'uint32');
    d_add = fread(fid_add, [3, params.chunk_size], 'uint32');
    if isempty(d_tr), break; end
    n_tr = size(d_tr, 2);

    b_tr = single(d_bs(1,:)')  ./ params.scale_factor;
    s_tr = max(single(d_add(1,:)'), 1e-6);

    id1 = d_tr(2,:) + 1;  id2 = d_tr(3,:) + 1;
    xs = single(geom.crystal_lut(id1, 1:3)'); xs = xs(:);
    xe = single(geom.crystal_lut(id2, 1:3)'); xe = xe(:);
    s1 = ones(n_tr, 1, 'single');

    z_low  = proj_forw_parallelproj_cuda(geom.libname, xs, xe, mu_low_img, ...
        geom.img_origin, geom.voxel_size, s1, n_tr, geom.img_dim);
    L_path = proj_forw_parallelproj_cuda(geom.libname, xs, xe, all_1_map, ...
        geom.img_origin, geom.voxel_size, s1, n_tr, geom.img_dim);

    z_low = min(max(z_low, 0), 50);
    u_i = b_tr .* exp(-z_low);
    w_i = max(u_i + s_tr, 1e-10);

    % Fisher 信息: mean(eta'^2) * sum(u^2/w * L^2)
    FI_tx = FI_tx + mean_eta_sq * sum(double(u_i.^2 ./ w_i .* L_path.^2));
end

fclose(fid_tr); fclose(fid_bs); fclose(fid_add);
end

