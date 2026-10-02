function image_recon = MLEM_parallelproj_TOF_bin_chunk_wAD_cuda(...
    libname, maxiter, image_recon, ...
    file_lm, file_td, ...
    sens, atten_map, ...                          % [新增] 传入当前的 atten_map
    PET_geom, NF_ax_data, NF_tr_data, ...         % [新增] 传入 Norm 需要的结构体
    crystal_lut, ...
    img_origin, voxel_size, img_dim, ...
    chunk_size, ...
    tofbin_width, sigma_tof, tofcenter_offset, n_sigmas_TOF, num_tof_bins, ... % [增加 num_tof_bins]
    lor_dependent_sigma_tof, lor_dependent_tofcenter_offset)

speed_of_light = 0.3; % mm/ps
epps = 1e-8;
img_dim = int32(img_dim);
n_voxels = prod(img_dim);

tof_bin_center_offset = int16(floor(double(num_tof_bins) / 2));  % = 37
tof_bin_max           = int16(num_tof_bins - 1);                 % = 74

for iter = 1 : maxiter
    disp(['[MLEM-Chunk] Processing iteration : ' num2str(iter)]);
    correction_accum = zeros(img_dim, 'single');
    correction_accum = correction_accum(:);

    fid_lm = fopen(file_lm, 'rb');
    fid_td = fopen(file_td, 'rb');
    if fid_lm == -1 || fid_td == -1
        error('无法打开 Listmode 或 TimeDiff 文件，请检查路径。');
    end

    chunk_counter = 1;
    while ~feof(fid_lm)
        % --- 1. 读取 Chunk ---
        d_lm = fread(fid_lm, [3, chunk_size], 'uint32');
        d_td = fread(fid_td, chunk_size, 'float');
        if isempty(d_lm) || isempty(d_td), break; end

        num_events = size(d_lm, 2);
        if length(d_td) ~= num_events
            num_events = length(d_td);
            d_lm = d_lm(:, 1:num_events);
            d_td = d_td(1:num_events);
        end

        % --- 2. 坐标与 TOF 计算 ---
        id1 = d_lm(2, :) + 1;
        id2 = d_lm(3, :) + 1;
        xstart = single(crystal_lut(id1, :)'); xstart = xstart(:);
        xend   = single(crystal_lut(id2, :)'); xend   = xend(:);

        % tof_bin_index = single(round(d_td .* (speed_of_light / (2 * tofbin_width)))); 
        tof_bin_index = int16(round(d_td .* (speed_of_light / (2 * tofbin_width)))) + tof_bin_center_offset;
        tof_bin_index = max(min(tof_bin_index, tof_bin_max), int16(0));

        sino_measured_chunk = ones(num_events, 1, 'single');
        % nlors_chunk = single(num_events);
        nlors_chunk = uint64(num_events);

        % --- 3. 物理校正因子计算 [新增] ---
        % (a) 获取 Attenuation 因子 a_i
        % sino_ones = ones(num_events, 1, 'single');
        sino_zeros = ones(num_events, 1, 'single');
        proj_mu = proj_forw_parallelproj_cuda(libname, xstart, xend, atten_map, img_origin, voxel_size, sino_zeros, num_events, img_dim);
        % 假设传入的 atten_map 已经是 mm^-1，如果不一致按需除以 10
        attn_factor = exp(-proj_mu);

        % (b) 获取 Normalization 因子 N_i
        LORs_frac = [d_lm(2, :)', d_lm(3, :)'];
        norm_factor = NF_generate_data(LORs_frac, crystal_lut, PET_geom, NF_ax_data, NF_tr_data);
        norm_inv = 1 ./ norm_factor;
        norm_inv(isnan(norm_inv))=0;
        norm_inv(isinf(norm_inv))=0;

        % (c) 获取 Additive 因子 r_it (Scatter & Randoms)
        [r_em, ~, ~] = get_additive_terms(id1, id2, num_events);

        % --- 4. TOF 正投影 (Forward Project) ---
        sino_fwd_act = proj_forw_parallelproj_TOF_bin_lm_cuda(...
            libname, xstart, xend, image_recon, img_origin, voxel_size, ...
            sino_measured_chunk, nlors_chunk, img_dim, ...
            tofbin_width, sigma_tof, tofcenter_offset, n_sigmas_TOF, ...
            tof_bin_index, num_tof_bins, lor_dependent_sigma_tof, lor_dependent_tofcenter_offset);

        % --- 5. 计算比值 (Ratio) [修改] ---
        % y_bar = a_i * N_i^(-1) * Px + r
        y_bar = attn_factor .* norm_inv .* sino_fwd_act + r_em;
        Ratio_factor_sino = sino_measured_chunk ./ (y_bar + epps);

        Ratio_factor_sino(isnan(Ratio_factor_sino)) = epps;
        Ratio_factor_sino(isinf(Ratio_factor_sino)) = epps;

        % --- 6. 反投影前赋予物理权重 [新增] ---
        % 因为反投影算子仅包含 P^T，我们需要反投影 (a_i * N_i^(-1) * Ratio)
        Ratio_to_backproject = Ratio_factor_sino .* attn_factor .* norm_inv;

        % --- 7. 反投影 (Back Project) ---
        back_chunk = proj_back_parallelproj_TOF_bin_lm_cuda(...
            libname, xstart, xend, image_recon, img_origin, voxel_size, ...
            Ratio_to_backproject, nlors_chunk, img_dim, ...
            tofbin_width, sigma_tof, tofcenter_offset, n_sigmas_TOF, ...
            tof_bin_index, num_tof_bins, lor_dependent_sigma_tof, lor_dependent_tofcenter_offset);

        correction_accum = correction_accum + back_chunk;
        chunk_counter = chunk_counter + 1;
    end
    fclose(fid_lm); fclose(fid_td);

    correction_accum(isnan(correction_accum)) = epps;
    correction_accum(isinf(correction_accum)) = epps;

    % --- 8. 图像更新 ---
    image_recon = image_recon ./ (sens + epps) .* correction_accum;
    image_recon(isnan(image_recon)) = epps;
    image_recon(isinf(image_recon)) = epps;
    image_recon(image_recon < epps) = epps;
end
end
