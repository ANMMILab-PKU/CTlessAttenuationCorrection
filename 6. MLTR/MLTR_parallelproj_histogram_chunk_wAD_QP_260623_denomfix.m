function atten_map_renew = MLTR_parallelproj_histogram_chunk_wAD_QP_260623_denomfix( ...
    libname, atten_map, filename_bs, filename_tr, filename_add, ...
    crystal_lut, img_dim, voxel_size, diemeter_scanner, alpha_step, chunk_size, scale_factor, beta_qp, factor)
% MLTR_PARALLELPROJ_HISTOGRAM_CHUNK_WAD_QP_260622_DENOMFIX
%
% Bug fix relative to ..._260622_fixed:
%   The denominator (surrogate curvature) weight was changed from c_i = u_i
%   to c_i = u_i / ybar_i. That broke the update:
%     - it is not a valid upper bound on the exact transmission curvature,
%     - it is dimensionally inconsistent (drops the "counts" dimension),
%       so the step is rescaled by ~ybar_i and the mu-map fails to evolve.
%   This version restores c_i = u_i (standard Nuyts MLTR curvature, a valid
%   upper bound because the exact curvature u_i*(1 - p_i*r_i/ybar_i^2) <= u_i).
%
% Statistical model:
%   p_i ~ Poisson( B_i * exp(-A_i*mu) + r_i )
%
% Update:
%   mu_j(new) = mu_j(old) + alpha *
%       [ sum_i L_ij u_i (1 - p_i/(u_i+r_i)) - beta_eff*dU/dmu_j ]
%       / [ D * sum_i L_ij u_i + beta_eff*d2U/dmu_j2 ]
%
% Scale convention (unchanged, matches 0523):
%   counts_bs_sub  = raw blank counts / scale_factor   (blank scaled DOWN to TX time)
%   counts_tr_sub  = raw TX all-events / prompt counts (NOT scaled)
%   counts_add_sub = raw TX additive  counts           (NOT scaled)

%% 1. Input defaults and type control
if nargin < 13 || isempty(beta_qp)
    beta_qp = single(0);
end
if nargin < 14 || isempty(factor)
    factor = single(1);
end

img_dim = int32(img_dim);
voxel_size = single(voxel_size);
atten_map = single(atten_map(:));

scale_factor = single(scale_factor);
beta_qp = single(beta_qp);
factor = single(factor);
alpha_step = single(alpha_step);
diemeter_scanner = single(diemeter_scanner);

if scale_factor <= 0
    error('scale_factor must be positive.');
end
if beta_qp < 0 || factor < 0
    error('beta_qp and factor must be non-negative.');
end

beta_eff = beta_qp .* factor;

n0 = single(img_dim(1));
n1 = single(img_dim(2));
n2 = single(img_dim(3));

img_origin = single([(-(n0 / 2 - 0.5)) * voxel_size(1), ...
                     (-(n1 / 2 - 0.5)) * voxel_size(2), ...
                     (-(n2 / 2 - 0.5)) * voxel_size(3)]);
clear n0 n1 n2

num_voxels = prod(double(img_dim));

Numerator_accum = zeros(num_voxels, 1, 'single');
Denominator_accum = zeros(num_voxels, 1, 'single');

zero_image = zeros(num_voxels, 1, 'single');

%% 2. Open histogram files
fid_bs  = fopen(filename_bs,  'rb');
fid_tr  = fopen(filename_tr,  'rb');
fid_add = fopen(filename_add, 'rb');

if fid_bs == -1 || fid_tr == -1 || fid_add == -1
    if fid_bs  ~= -1, fclose(fid_bs);  end
    if fid_tr  ~= -1, fclose(fid_tr);  end
    if fid_add ~= -1, fclose(fid_add); end
    error('Error: Could not open histogram data files.');
end

disp(['    [MLTR-QP] Stream processing started. beta_qp = ' ...
      num2str(double(beta_qp)) ', factor = ' num2str(double(factor)) ...
      ', beta_eff = ' num2str(double(beta_eff)) ...
      ', scale_factor = ' num2str(double(scale_factor))]);

%% 3. Chunk-wise MLTR accumulation
try
    while ~feof(fid_tr) && ~feof(fid_bs) && ~feof(fid_add)

        % --- 3.1 Read one data chunk from each unified histogram ---
        data_bs_chunk  = fread(fid_bs,  [3, chunk_size], 'uint32');
        data_tr_chunk  = fread(fid_tr,  [3, chunk_size], 'uint32');
        data_add_chunk = fread(fid_add, [3, chunk_size], 'uint32');

        if isempty(data_tr_chunk) || isempty(data_bs_chunk) || isempty(data_add_chunk)
            break;
        end

        num_events = min([size(data_bs_chunk, 2), ...
                          size(data_tr_chunk, 2), ...
                          size(data_add_chunk, 2)]);

        if num_events <= 0
            break;
        end

        data_bs_chunk  = data_bs_chunk(:,  1:num_events);
        data_tr_chunk  = data_tr_chunk(:,  1:num_events);
        data_add_chunk = data_add_chunk(:, 1:num_events);

        % --- 3.2 Time-scale alignment (blank scaled down to TX time) ---
        counts_bs_sub  = single(data_bs_chunk(1, :)') ./ scale_factor;
        counts_tr_sub  = single(data_tr_chunk(1, :)');
        counts_add_sub = single(data_add_chunk(1, :)');

        id1 = double(data_tr_chunk(2, :)) + 1;
        id2 = double(data_tr_chunk(3, :)) + 1;

        xstart_sub = single(crystal_lut(id1, :)');
        xstart_sub = xstart_sub(:);

        xend_sub = single(crystal_lut(id2, :)');
        xend_sub = xend_sub(:);

        sino_init_sub = ones(num_events, 1, 'single');

        % --- 3.3 Forward projection z = sum_j L_ij * mu_j ---
        sino_fwd_sub = proj_forw_parallelproj( ...
            libname, xstart_sub, xend_sub, atten_map, ...
            img_origin, voxel_size, sino_init_sub, num_events, img_dim);

        % --- 3.4 Additive transmission model ---
        u_expected = counts_bs_sub .* exp((-1) .* sino_fwd_sub);
        y_bar = u_expected + counts_add_sub;
        y_bar_safe = max(y_bar, single(1e-9));

        % --- 3.5 MLTR numerator and denominator terms ---
        % Numerator data term:  u_i * (1 - p_i / (u_i + r_i))
        term_gradient = u_expected .* (1 - counts_tr_sub ./ y_bar_safe);

        % Denominator (surrogate curvature) term:  c_i = u_i
        % Valid upper bound on exact curvature u_i*(1 - p_i*r_i/ybar_i^2) <= u_i.
        % (FIX: previous version used u_i ./ y_bar_safe, which is wrong.)
        term_weight = u_expected;

        % --- 3.6 Backproject and accumulate ---
        Numerator_update = proj_back_parallelproj( ...
            libname, xstart_sub, xend_sub, zero_image, ...
            img_origin, voxel_size, term_gradient, num_events, img_dim);
        Numerator_accum = Numerator_accum + single(Numerator_update(:));

        Denominator_update = proj_back_parallelproj( ...
            libname, xstart_sub, xend_sub, zero_image, ...
            img_origin, voxel_size, term_weight, num_events, img_dim);
        Denominator_accum = Denominator_accum + single(Denominator_update(:));
    end

    fclose(fid_bs);
    fclose(fid_tr);
    fclose(fid_add);

catch ME
    if fid_bs  ~= -1, fclose(fid_bs);  end
    if fid_tr  ~= -1, fclose(fid_tr);  end
    if fid_add ~= -1, fclose(fid_add); end
    rethrow(ME);
end

%% 4. Add QP regularization and update mu-map
disp('    [MLTR-QP] Applying update...');

Numerator_final = Numerator_accum;
Denominator_final = diemeter_scanner .* Denominator_accum;

if beta_eff > 0
    atten_img = reshape(atten_map, double(img_dim(:)'));
    [reg_grad, reg_curv] = compute_quadratic_prior_260622_fixed(atten_img, beta_eff);
    Numerator_final = Numerator_final - reg_grad(:);
    Denominator_final = Denominator_final + reg_curv(:);
end

Denominator_final = Denominator_final + single(1e-9);

update_miu = alpha_step .* (Numerator_final ./ Denominator_final);
update_miu(isnan(update_miu) | isinf(update_miu)) = 0;

atten_map_renew = atten_map + update_miu;

% Non-negative constraint for mu-map.
atten_map_renew(atten_map_renew < 0) = 0;
atten_map_renew = single(atten_map_renew);

end
