function atten_map_renew = ML_attn_sps_parallelproj_lm_norm_chunk_cuda( ...
    libname, atten_map, activity_map, ...
    file_lm, ...
    crystal_lut, PET_geom, NF_ax_data, NF_tr_data, ...
    img_origin, voxel_size, img_dim, chunk_size, ...
    num_parts, ring_num, crystal_per_ring, ...
    min_crystal_difference, max_ring_difference)
% Attenuation Maximum Likelihood using Symmetric-Pixel-Search
%
% Chunk + CUDA + unified Norm-data version.
%
% Main changes:
%   1) The listmode contribution nor_yi is accumulated chunk-by-chunk from
%      file_lm, following the same fread([3, chunk_size], 'uint32') pattern
%      used in MLEM_parallelproj_TOF_bin_chunk_wAD_cuda.
%   2) Forward/back projectors use CUDA wrappers:
%          proj_forw_parallelproj_cuda
%          proj_back_parallelproj_cuda
%   3) Normalization uses loaded NF_ax_data / NF_tr_data through
%      NF_generate_data, matching the MLEM interface style.
%
% Required listmode file format:
%   file_lm stores uint32 records readable as [3, chunk_size].
%   Crystal IDs are read the same way as in the reference MLEM:
%       id1 = d_lm(2,:) + 1;
%       id2 = d_lm(3,:) + 1;
%
% Unit note:
%   The original attenuation convention is preserved:
%       attenuation_factor = exp(-sino_fwd ./ 10)
%   If atten_map is already in mm^-1 and the projector line integral is
%   already in matching mm units, change this to exp(-sino_fwd).

epps = single(1e-8);

img_dim = int32(img_dim);
chunk_size = double(chunk_size);
n_voxels = prod(double(img_dim));

atten_map = single(atten_map(:));
activity_map = single(activity_map(:));
atten_map_renew = atten_map;

% CUDA back projector only needs this as an output buffer template.
image_recon = zeros(n_voxels, 1, 'single');

%% ------------------------------------------------------------------------
% 1. gradient_yi: listmode measured data backprojection, accumulated by chunk
% -------------------------------------------------------------------------

disp('[ML-attn-SPS-Chunk-CUDA] Accumulating nor_yi from listmode chunks ...');

nor_yi = zeros(n_voxels, 1, 'single');

fid_lm = fopen(file_lm, 'rb');
if fid_lm == -1
    error('无法打开 Listmode 文件，请检查路径: %s', file_lm);
end
cleanup_lm = onCleanup(@() fclose(fid_lm));

chunk_counter = 1;
while ~feof(fid_lm)
    d_lm = fread(fid_lm, [3, chunk_size], 'uint32');
    if isempty(d_lm)
        break;
    end

    num_events = size(d_lm, 2);
    if num_events == 0
        break;
    end

    id1 = double(d_lm(2, :)) + 1;
    id2 = double(d_lm(3, :)) + 1;

    xstart = single(crystal_lut(id1, :)');
    xend   = single(crystal_lut(id2, :)');
    xstart = xstart(:);
    xend   = xend(:);

    sino_yi_projection = ones(num_events, 1, 'single');
    nlors_lmf = uint64(num_events);

    nor_yi_chunk = proj_back_parallelproj_cuda( ...
        libname, xstart, xend, image_recon, ...
        img_origin, voxel_size, sino_yi_projection, nlors_lmf, img_dim);

    nor_yi = nor_yi + single(nor_yi_chunk(:));

    clear d_lm id1 id2 xstart xend sino_yi_projection nor_yi_chunk

    if mod(chunk_counter, 20) == 0
        disp(['[ML-attn-SPS-Chunk-CUDA] Finished LM chunk ' num2str(chunk_counter)]);
    end
    chunk_counter = chunk_counter + 1;
end

clear cleanup_lm fid_lm

%% ------------------------------------------------------------------------
% 2. Expected data terms over all valid LORs, split by num_parts
% -------------------------------------------------------------------------

all_1_map = ones(n_voxels, 1, 'single');

denr_all    = zeros(n_voxels, 1, 'single');
nor_est_all = zeros(n_voxels, 1, 'single');

for i_part = 1 : num_parts

    % disp(['[ML-attn-SPS-Chunk-CUDA] Calculate valid LOR part ' num2str(i_part) ' / ' num2str(num_parts)]);

    sino_id = Get_valid_lor_calculate_sens_single_layer( ...
        ring_num, crystal_per_ring, ...
        min_crystal_difference, max_ring_difference, ...
        num_parts, i_part);

    if isempty(sino_id)
        continue;
    end

    coor_sinoid1 = crystal_lut(sino_id(:, 1), :);
    coor_sinoid2 = crystal_lut(sino_id(:, 2), :);

    xstart = single(coor_sinoid1');
    xend   = single(coor_sinoid2');
    xstart = xstart(:);
    xend   = xend(:);

    clear coor_sinoid1 coor_sinoid2

    nlors_miumap = uint64(size(sino_id, 1));
    sino_init = ones(size(sino_id, 1), 1, 'single');

    % attenuation factor
    sino_fwd = proj_forw_parallelproj_cuda( ...
        libname, xstart, xend, atten_map_renew, ...
        img_origin, voxel_size, sino_init, nlors_miumap, img_dim);

    attenuation_factor = exp((-1) * sino_fwd );

    % normalization factor, using the same loaded-data style as MLEM
    LORs_frac = sino_id - 1;
    normalization_factor = NF_generate_data( ...
        LORs_frac, crystal_lut, PET_geom, NF_ax_data, NF_tr_data);
    normalization_factor = single(normalization_factor(:));

    norm_inv = 1 ./ normalization_factor;
    norm_inv(isnan(norm_inv) | isinf(norm_inv)) = 0;

    % activity projection
    projection_activityMap = proj_forw_parallelproj_cuda( ...
        libname, xstart, xend, activity_map, ...
        img_origin, voxel_size, sino_init, nlors_miumap, img_dim);

    estimation = attenuation_factor .* projection_activityMap .* norm_inv;
    estimation(isnan(estimation) | isinf(estimation)) = 0;

    % gradient_est
    nor_est = proj_back_parallelproj_cuda( ...
        libname, xstart, xend, image_recon, ...
        img_origin, voxel_size, estimation, nlors_miumap, img_dim);

    % sum_G_trans
    sum_G_trans = proj_forw_parallelproj_cuda( ...
        libname, xstart, xend, all_1_map, ...
        img_origin, voxel_size, sino_init, nlors_miumap, img_dim);

    % est_sum_G_trans
    est_sum_G_trans = estimation .* sum_G_trans;
    est_sum_G_trans(isnan(est_sum_G_trans) | isinf(est_sum_G_trans)) = 0;

    % curvature
    denr = proj_back_parallelproj_cuda( ...
        libname, xstart, xend, image_recon, ...
        img_origin, voxel_size, est_sum_G_trans, nlors_miumap, img_dim);

    nor_est_all = nor_est_all + single(nor_est(:));
    denr_all    = denr_all    + single(denr(:));

    clear sino_id LORs_frac xstart xend sino_init
    clear sino_fwd attenuation_factor normalization_factor norm_inv
    clear projection_activityMap estimation nor_est
    clear sum_G_trans est_sum_G_trans denr
end

%% ------------------------------------------------------------------------
% 3. SPS update
% -------------------------------------------------------------------------

nor_all = nor_est_all - nor_yi;

update_miu = nor_all ./ (denr_all + epps);
update_miu(isnan(update_miu) | isinf(update_miu)) = 0;

atten_map_renew = atten_map_renew + update_miu;

% miu-map non-negative
atten_map_renew(atten_map_renew < 0) = 0;

end
