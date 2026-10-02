function atten_map_renew = MLAA_TX_parallelproj_histo_chunk_data_wAD_LCcul_3MEP( ...
    libname, atten_map, activity_map, ...
    file_em, ...
    file_tr_202, file_bs_202, file_tr_202_add, ... 
    file_tr_307, file_bs_307, file_tr_307_add, ... 
    PET_lut, crystal_lut, img_dim, voxel_size, ...
    alpha_EM_511, alpha_TX_202, alpha_TX_307, ...
    beta_reg, ...
    chunk_size, ...
    tx_mode, ...
    PET_geom, NF_ax_data, NF_tr_data, ...
    num_parts, ring_num, crystal_per_ring, ...
    min_crystal_difference, max_ring_difference, scale_factor)

% MLAA_TX_PARALLELPROJ_HISTO_CHUNK (Multi-Energy Mapping Version)
% 整合分段线性映射模型，解决骨骼高估问题

%% ── 1. Shared setup ───────────────────────────────────────────────────────

img_dim    = int32(img_dim);
n0 = single(img_dim(1)); n1 = single(img_dim(2)); n2 = single(img_dim(3));
img_origin = single([ (-(n0/2-0.5))*voxel_size(1), (-(n1/2-0.5))*voxel_size(2), (-(n2/2-0.5))*voxel_size(3) ]);

use_202 = (tx_mode == 1 || tx_mode == 3);
use_307 = (tx_mode == 2 || tx_mode == 3);

Nvox       = prod(img_dim);
zero_image = single(zeros(Nvox, 1));
all_1_map  = single(ones( Nvox, 1));

% scale_factor = (10000 * 10) / (240 * 10); % Blank/Phantom time scale	% ≈41.67
% scale_factor =  (5000 * 7) / (600 * 2); % ≈ 29.17



% --- [关键修改]：执行低能映射与导数计算 ---
if use_202
    [mu_202_map, deriv_202_map] = map_mu511_to_low_energy_3segLM(atten_map, 202);
end
if use_307
    [mu_307_map, deriv_307_map] = map_mu511_to_low_energy_3segLM(atten_map, 307);
end

s_mu_em    = zeros(Nvox, 1, 'single');
c_em_total = zeros(Nvox, 1, 'single');

%% ── 2. Pass 0: Emission Sensitivity (r=0 for True Data) ──────────────────

if alpha_EM_511 > 0
    fprintf('  Pass 0: sensitivity loop...\n');
    for i_part = 1 : num_parts
        sino_id = Get_valid_lor_calculate_sens_single_layer(ring_num, crystal_per_ring, min_crystal_difference, max_ring_difference, num_parts, i_part);
        xs_p = single(crystal_lut(sino_id(:,1), :)'); xs_p = xs_p(:);
        xe_p = single(crystal_lut(sino_id(:,2), :)'); xe_p = xe_p(:);
        nlors_p = size(sino_id, 1);
        s1 = single(ones(nlors_p, 1));

        proj_mu_p  = proj_forw_parallelproj_cuda(libname, xs_p, xe_p, atten_map,    img_origin, voxel_size, s1, nlors_p, img_dim);
        proj_act_p = proj_forw_parallelproj_cuda(libname, xs_p, xe_p, activity_map, img_origin, voxel_size, s1, nlors_p, img_dim);
        L_path_p   = proj_forw_parallelproj_cuda(libname, xs_p, xe_p, all_1_map,    img_origin, voxel_size, s1, nlors_p, img_dim);

        nf_p = NF_generate_data(sino_id-1, crystal_lut, PET_geom, NF_ax_data, NF_tr_data);
        phi_p = proj_act_p ./ nf_p;
        phi_p(isnan(phi_p) | isinf(phi_p)) = 0;
        
        r_p = zeros(nlors_p, 1, 'single'); % 预校正数据 r=0
        z_p = proj_mu_p;
        y_bar_p = phi_p .* exp(-z_p) + r_p;
        
        C_bg_p = calc_optimal_curvature_double_opt_0p001(z_p, zeros(nlors_p,1,'single'), y_bar_p, r_p);
        psi_p = phi_p .* exp(-z_p);

        s_mu_em    = s_mu_em    + proj_back_parallelproj_cuda(libname, xs_p, xe_p, zero_image, img_origin, voxel_size, psi_p, nlors_p, img_dim);
        c_em_total = c_em_total + proj_back_parallelproj_cuda(libname, xs_p, xe_p, zero_image, img_origin, voxel_size, C_bg_p .* L_path_p, nlors_p, img_dim); 
    end
end

Numerator_accum   = alpha_EM_511 .* s_mu_em;
Denominator_accum = alpha_EM_511 .* c_em_total;

%% ── 3. Loop A: Emission Term B (Sparse y_em > 0) ─────────────────────────

if alpha_EM_511 > 0
    fprintf('  Loop A: emission Term B...\n');
    fid_em = fopen(file_em, 'rb');
    while ~feof(fid_em)
        d_em = fread(fid_em, [3, chunk_size], 'uint32');
        if isempty(d_em), break; end
        n_em = size(d_em, 2);
        xs = single(crystal_lut(d_em(2,:)+1, :)'); xs = xs(:);
        xe = single(crystal_lut(d_em(3,:)+1, :)'); xe = xe(:);
        y_em = single(d_em(1,:)');
        s1 = single(ones(n_em, 1));

        proj_mu = proj_forw_parallelproj_cuda(libname, xs, xe, atten_map, img_origin, voxel_size, s1, n_em, img_dim);
        proj_act = proj_forw_parallelproj_cuda(libname, xs, xe, activity_map, img_origin, voxel_size, s1, n_em, img_dim);
        nf_em = NF_generate_data([d_em(2,:)', d_em(3,:)'], crystal_lut, PET_geom, NF_ax_data, NF_tr_data);

        phi_em = proj_act ./ nf_em;
        phi_em(isnan(phi_em) | isinf(phi_em)) = 0;
        z_em = proj_mu;
        y_bar_em = phi_em .* exp(-z_em); 

        data_weight = (phi_em .* exp(-z_em)) .* (y_em ./ (y_bar_em + 1e-9));
        Numerator_accum = Numerator_accum - alpha_EM_511 .* proj_back_parallelproj_cuda(libname, xs, xe, zero_image, img_origin, voxel_size, data_weight, n_em, img_dim);

        L_path_em = proj_forw_parallelproj_cuda(libname, xs, xe, all_1_map, img_origin, voxel_size, s1, n_em, img_dim);
        C_meas_em = calc_optimal_curvature_double_opt_0p001(z_em, y_em, y_bar_em, zeros(n_em,1,'single'));
        C_bg_em   = calc_optimal_curvature_double_opt_0p001(z_em, zeros(n_em,1,'single'), y_bar_em, zeros(n_em,1,'single'));
        Denominator_accum = Denominator_accum + alpha_EM_511 .* proj_back_parallelproj_cuda(libname, xs, xe, zero_image, img_origin, voxel_size, (C_meas_em - C_bg_em) .* L_path_em, n_em, img_dim);
    end
    fclose(fid_em);
end

%% ── 4. Loop B: TX 202 keV (Bilinear Mapping) ──────────────────

if use_202
    fprintf('  Loop B: TX-202 Unified (Mapping)...\n');
    fid_tr = fopen(file_tr_202, 'rb'); fid_bs = fopen(file_bs_202, 'rb'); fid_add = fopen(file_tr_202_add, 'rb');
    while ~feof(fid_tr)
        d_tr = fread(fid_tr, [3, chunk_size], 'uint32'); d_bs = fread(fid_bs, [3, chunk_size], 'uint32'); d_add = fread(fid_add, [3, chunk_size], 'uint32');
        if isempty(d_tr), break; end
        n_tr = size(d_tr, 2);
        xs = single(crystal_lut(d_tr(2,:)+1, :)'); xs = xs(:);
        xe = single(crystal_lut(d_tr(3,:)+1, :)'); xe = xe(:);
        y_tr = single(d_tr(1,:)'); b_tr = single(d_bs(1,:)') ./ scale_factor; s_tr = single(d_add(1,:)');

        % 使用映射后的低能衰减图进行前向投影
        z_202 = proj_forw_parallelproj_cuda(libname, xs, xe, mu_202_map, img_origin, voxel_size, single(ones(n_tr,1)), n_tr, img_dim);
        y_signal = b_tr .* exp(-z_202);
        y_bar = y_signal + s_tr;

        % LOR 空间梯度与曲率
        grad_lor = alpha_TX_202 .* (1 - y_tr ./ (y_bar + 1e-9)) .* y_signal;
        curv_lor = alpha_TX_202 .* calc_optimal_curvature_double_opt_0p001(z_202, y_tr, y_bar, s_tr) .* proj_forw_parallelproj_cuda(libname, xs, xe, all_1_map, img_origin, voxel_size, single(ones(n_tr,1)), n_tr, img_dim);

        % 链式法则：乘以体素级导数 deriv_eta
        Numerator_accum   = Numerator_accum   + deriv_202_map .* proj_back_parallelproj_cuda(libname, xs, xe, zero_image, img_origin, voxel_size, grad_lor, n_tr, img_dim);
        Denominator_accum = Denominator_accum + (deriv_202_map.^2) .* proj_back_parallelproj_cuda(libname, xs, xe, zero_image, img_origin, voxel_size, curv_lor, n_tr, img_dim);
    end
    fclose(fid_tr); fclose(fid_bs); fclose(fid_add);
end

%% ── 5. Loop C: TX 307 keV (Bilinear Mapping) ──────────────────

if use_307
    fprintf('  Loop C: TX-307 Unified (Mapping)...\n');
    fid_tr = fopen(file_tr_307, 'rb'); fid_bs = fopen(file_bs_307, 'rb'); fid_add = fopen(file_tr_307_add, 'rb');
    while ~feof(fid_tr)
        d_tr = fread(fid_tr, [3, chunk_size], 'uint32'); d_bs = fread(fid_bs, [3, chunk_size], 'uint32'); d_add = fread(fid_add, [3, chunk_size], 'uint32');
        if isempty(d_tr), break; end
        n_tr = size(d_tr, 2);
        xs = single(crystal_lut(d_tr(2,:)+1, :)'); xs = xs(:);
        xe = single(crystal_lut(d_tr(3,:)+1, :)'); xe = xe(:);
        y_tr = single(d_tr(1,:)'); b_tr = single(d_bs(1,:)') ./ scale_factor; s_tr = single(d_add(1,:)');

        z_307 = proj_forw_parallelproj_cuda(libname, xs, xe, mu_307_map, img_origin, voxel_size, single(ones(n_tr,1)), n_tr, img_dim);
        y_signal = b_tr .* exp(-z_307);
        y_bar = y_signal + s_tr;

        grad_lor = alpha_TX_307 .* (1 - y_tr ./ (y_bar + 1e-9)) .* y_signal;
        curv_lor = alpha_TX_307 .* calc_optimal_curvature_double_opt_0p001(z_307, y_tr, y_bar, s_tr) .* proj_forw_parallelproj_cuda(libname, xs, xe, all_1_map, img_origin, voxel_size, single(ones(n_tr,1)), n_tr, img_dim);

        Numerator_accum   = Numerator_accum   + deriv_307_map .* proj_back_parallelproj_cuda(libname, xs, xe, zero_image, img_origin, voxel_size, grad_lor, n_tr, img_dim);
        Denominator_accum = Denominator_accum + (deriv_307_map.^2) .* proj_back_parallelproj_cuda(libname, xs, xe, zero_image, img_origin, voxel_size, curv_lor, n_tr, img_dim);
    end
    fclose(fid_tr); fclose(fid_bs); fclose(fid_add);
end

%% ── 6. Final Update ──────────────────────────────────────────────────────

if beta_reg > 0
    [reg_grad, reg_curv] = calc_regularization_quadratic(atten_map, beta_reg, img_dim, voxel_size);
    Numerator_accum = Numerator_accum - reg_grad;
    Denominator_accum = Denominator_accum + reg_curv;
end

update_step = Numerator_accum ./ (Denominator_accum + 1e-9);
atten_map_renew = max(atten_map(:) + update_step, 0);

end
