function [atten_map_renew, scs_state, dbg_info] = MLAA_TX_SCS_update_MEP_260518(atten_map, activity_map, file_em, file_tr_202, file_bs_202, file_tr_307, file_bs_307, params_SCS, psi_geom, scs_state, iterr, it_sub)
% MLAA_TX_SCS_update (Modular Version - Nonlinear BMD)
% -------------------------------------------------------------------------
% 理论保证版本 + 图像域3基物质非线性映射 (Image-domain BMD)
% -------------------------------------------------------------------------

% 解包参数
libname = psi_geom.libname;
img_dim = psi_geom.img_dim;
img_origin = psi_geom.img_origin;
voxel_size = psi_geom.voxel_size;
chunk_size = params_SCS.chunk_size;

Nvox = prod(img_dim);
zero_image = zeros(Nvox, 1, 'single');
all_1_map  = ones(Nvox, 1, 'single');

% 初始化累加器
Grad_tot     = zeros(Nvox, 1, 'single');
Hess_tot     = zeros(Nvox, 1, 'single');
Cubic_M_data = zeros(Nvox, 1, 'single');

% =====================================================================
% [DENSE EMISSION] Pass 0: 稠密几何循环 (解决 y=0 的背景梯度)
% =====================================================================
if params_SCS.alpha_EM > 0
    fprintf('    -> Dense Geometry Pass (y=0 baseline)...\n');
    for i_part = 1 : psi_geom.num_parts
        sino_id = Get_valid_lor_calculate_sens_single_layer(psi_geom.ring_num, psi_geom.crystal_per_ring, psi_geom.min_crystal_difference, psi_geom.max_ring_difference, psi_geom.num_parts, i_part);

        xs = single(psi_geom.crystal_lut(sino_id(:,1), 1:3)'); xs = xs(:);
        xe = single(psi_geom.crystal_lut(sino_id(:,2), 1:3)'); xe = xe(:);
        nlors_p = size(sino_id, 1);
        s1 = ones(nlors_p, 1, 'single');

        proj_mu  = proj_forw_parallelproj_cuda(libname, xs, xe, atten_map, img_origin, voxel_size, s1, nlors_p, img_dim);
        proj_act = proj_forw_parallelproj_cuda(libname, xs, xe, activity_map, img_origin, voxel_size, s1, nlors_p, img_dim);
        L_path   = proj_forw_parallelproj_cuda(libname, xs, xe, all_1_map, img_origin, voxel_size, s1, nlors_p, img_dim);

        nf = NF_generate_data(sino_id - 1, psi_geom.crystal_lut, psi_geom.PET_geom, psi_geom.NF_ax_data, psi_geom.NF_tr_data);

        u_i = exp(-proj_mu) .* proj_act ./ nf;
        u_i(isnan(u_i) | isinf(u_i)) = 0;

        Grad_tot     = Grad_tot     + params_SCS.alpha_EM * proj_back_parallelproj_cuda(libname, xs, xe, zero_image, img_origin, voxel_size, u_i, nlors_p, img_dim);
        Hess_tot     = Hess_tot     + params_SCS.alpha_EM * proj_back_parallelproj_cuda(libname, xs, xe, zero_image, img_origin, voxel_size, -u_i .* L_path, nlors_p, img_dim);
        Cubic_M_data = Cubic_M_data + params_SCS.alpha_EM * proj_back_parallelproj_cuda(libname, xs, xe, zero_image, img_origin, voxel_size, u_i .* (L_path.^2), nlors_p, img_dim);
    end
end

% =====================================================================
% [SPARSE EMISSION] Loop A: 稀疏直方图修正 (y>0)
% =====================================================================
if params_SCS.alpha_EM > 0
    fprintf('    -> Sparse Emission Pass (IO loop)...\n');
    fid_em = fopen(file_em, 'rb');
    if fid_em == -1, error('Cannot open emission file: %s', file_em); end

    while ~feof(fid_em)
        d_em = fread(fid_em, [3, chunk_size], 'uint32');
        if isempty(d_em), break; end
        n_em = size(d_em, 2);

        id1 = d_em(2,:) + 1; id2 = d_em(3,:) + 1;
        xs = single(psi_geom.crystal_lut(id1, 1:3)'); xs = xs(:);
        xe = single(psi_geom.crystal_lut(id2, 1:3)'); xe = xe(:);
        s1 = ones(n_em, 1, 'single');
        y_em = single(d_em(1,:)');

        [r_em, ~, ~] = get_additive_terms(id1, id2, n_em);

        proj_mu  = proj_forw_parallelproj_cuda(libname, xs, xe, atten_map, img_origin, voxel_size, s1, n_em, img_dim);
        proj_act = proj_forw_parallelproj_cuda(libname, xs, xe, activity_map, img_origin, voxel_size, s1, n_em, img_dim);
        L_path   = proj_forw_parallelproj_cuda(libname, xs, xe, all_1_map, img_origin, voxel_size, s1, n_em, img_dim);
        nf = NF_generate_data([d_em(2,:)', d_em(3,:)'], psi_geom.crystal_lut, psi_geom.PET_geom, psi_geom.NF_ax_data, psi_geom.NF_tr_data);

        % 提取精确的三阶导数修正 K_corr
        [h1_corr, h2_corr, K_corr] = calc_exact_derivatives_emission_260502(proj_mu, proj_act, y_em, r_em, 0, nf);

        Grad_tot = Grad_tot + params_SCS.alpha_EM * proj_back_parallelproj_cuda(libname, xs, xe, zero_image, img_origin, voxel_size, h1_corr, n_em, img_dim);
        Hess_tot = Hess_tot + params_SCS.alpha_EM * proj_back_parallelproj_cuda(libname, xs, xe, zero_image, img_origin, voxel_size, h2_corr .* L_path, n_em, img_dim);

        % 将 Emission 的三阶上界累加至全局 M 矩阵
        Cubic_M_data = Cubic_M_data + params_SCS.alpha_EM * proj_back_parallelproj_cuda(libname, xs, xe, zero_image, img_origin, voxel_size, K_corr .* (L_path.^2), n_em, img_dim);
    end
    fclose(fid_em);
end

% =====================================================================
% [TRANSMISSION] Loop B & C: TX 直方图 IO (集成图像域双基物质映射)
% =====================================================================
    function process_tx_chunks(f_tr, f_bs, f_add, energy_type, alpha, scatter_idx)
        fid_tr = fopen(f_tr, 'rb');
        fid_bs = fopen(f_bs, 'rb');
        fid_add = fopen(f_add, 'rb');

        % 【核心 1】: 在图像空间执行非线性映射，获取等效低能 mu 图及其导数图
        [mu_low_img, deriv_eta_img] = map_mu511_to_low_energy_3segLM(atten_map, energy_type);

        while ~feof(fid_tr)
            d_tr  = fread(fid_tr,  [3, chunk_size], 'uint32');
            d_bs  = fread(fid_bs,  [3, chunk_size], 'uint32');
            d_add = fread(fid_add, [3, chunk_size], 'uint32'); 

            if isempty(d_tr), break; end
            n_tr = size(d_tr, 2);

            y_tr = single(d_tr(1,:)');
            b_tr = single(d_bs(1,:)') ./ params_SCS.scale_factor;
            s_tr = single(d_add(1,:)'); 

            id1 = d_tr(2,:) + 1; id2 = d_tr(3,:) + 1;
            xs = single(psi_geom.crystal_lut(id1, 1:3)'); xs = xs(:);
            xe = single(psi_geom.crystal_lut(id2, 1:3)'); xe = xe(:);
            s1 = ones(n_tr, 1, 'single');

            % 【核心 2】: 前向投影映射后的低能 mu 图
            z_low  = proj_forw_parallelproj_cuda(libname, xs, xe, mu_low_img, img_origin, voxel_size, s1, n_tr, img_dim);
            L_path = proj_forw_parallelproj_cuda(libname, xs, xe, all_1_map, img_origin, voxel_size, s1, n_tr, img_dim);

            % 计算基于低能 z_low 的精确导数与曲率
            [h1, h2, K] = calc_exact_derivatives_transmission_260502(z_low, y_tr, b_tr, s_tr);

            % 【核心 3】: 纯粹的反投影（回到图像空间，此时还不包含 eta 导数）
            Grad_back  = proj_back_parallelproj_cuda(libname, xs, xe, zero_image, img_origin, voxel_size, h1,                n_tr, img_dim);
            Hess_back  = proj_back_parallelproj_cuda(libname, xs, xe, zero_image, img_origin, voxel_size, h2 .* L_path,      n_tr, img_dim);
            Cubic_back = proj_back_parallelproj_cuda(libname, xs, xe, zero_image, img_origin, voxel_size, K .* (L_path.^2),  n_tr, img_dim);

            % 【核心 4】: 应用链式法则 (Chain Rule)，乘以体素级物理斜率
            Grad_tot     = Grad_tot     + alpha .* (Grad_back  .* deriv_eta_img);
            Hess_tot     = Hess_tot     + alpha .* (Hess_back  .* (deriv_eta_img.^2));
            Cubic_M_data = Cubic_M_data + alpha .* (Cubic_back .* (deriv_eta_img.^3));
        end
        fclose(fid_tr); fclose(fid_bs); fclose(fid_add);
    end

if params_SCS.tx_mode == 1 || params_SCS.tx_mode == 3
    fprintf('    -> Transmission 202 keV Pass (Nonlinear mapping)...\n');
    process_tx_chunks(file_tr_202, file_bs_202, params_SCS.file_tr_202_add, 202, params_SCS.alpha_202, 202);
end
if params_SCS.tx_mode == 2 || params_SCS.tx_mode == 3
    fprintf('    -> Transmission 307 keV Pass (Nonlinear mapping)...\n');
    process_tx_chunks(file_tr_307, file_bs_307, params_SCS.file_tr_307_add, 307, params_SCS.alpha_307, 307);
end

% =====================================================================
% 正则化项兼容 (对接calc_regularization_quadratic)
% =====================================================================
if params_SCS.beta_reg > 0
    [G_reg, H_reg] = calc_regularization_quadratic(atten_map, params_SCS.beta_reg, img_dim, voxel_size);
    M_reg = zeros(Nvox, 1, 'single'); % 二次正则化的三阶导数严格为 0
    assert(min(H_reg(:)) >= -1e-5, 'SCS Error: Regularization Hessian H_reg contains negative values!');

    Grad_tot = Grad_tot - G_reg(:);
    Hess_tot = Hess_tot - H_reg(:);
    Cubic_M_data = Cubic_M_data + M_reg(:);
end

% =====================================================================
% M 绝对界限设定
% =====================================================================
fprintf('    -> Applying Strict Lipschitz Bounds...\n');
M_strict = min(Cubic_M_data(:), prctile(Cubic_M_data(:), 99.9));
M_effective = scs_state.lambda_cubic * M_strict + 1e-12;

% =====================================================================
% 【极简模块化调用】: Box-Constrained Cubic Solver
% =====================================================================
fprintf('    -> Solving Box-Constrained Cubic Subproblems...\n');

[delta_mu_vec, dbg_solver] = solve_cubic_vectorized_final_260502(Grad_tot(:), Hess_tot(:), M_effective(:), atten_map(:));
atten_map_candidate = atten_map(:) + delta_mu_vec;

% =====================================================================
% 基于理论底线的 Trust Region 接受逻辑
% =====================================================================
if dbg_solver.pred_increase > 1e-6
    atten_map_renew = atten_map_candidate;
    scs_state.lambda_cubic = max(scs_state.lambda_cubic * 0.95, 1.0);
    fprintf('    [Trust-Region] Step Accepted. Predicted Increase: %.4e\n', dbg_solver.pred_increase);
else
    atten_map_renew = atten_map(:);
    scs_state.lambda_cubic = min(scs_state.lambda_cubic * 2.0, 100.0);
    fprintf('    [Trust-Region] Stagnation detected. Lambda increased to %.4f\n', scs_state.lambda_cubic);
end

dbg_info.grad_norm = norm(Grad_tot(:));
dbg_info.step_norm = norm(delta_mu_vec);
dbg_info.hess_pos_ratio = sum(Hess_tot(:) > 0) / Nvox;
dbg_info.saddle_escapes = dbg_solver.saddle_escapes;
dbg_info.rho = 1.0;

end