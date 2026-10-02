clc
clear

initial_load = 0;

disp(mfilename);

%% 1. Reconstruction Parameters

maxit = 50;
sub_it = 4;
it_local = 1;

chunk_size = 10000000; % 10M events per chunk

%% 2. Image Geometry

img_dim = int32([350, 350, 138]);
voxel_size = single([2, 2, 2]);

n0 = single(img_dim(1));
n1 = single(img_dim(2));
n2 = single(img_dim(3));
img_origin = single([ ...
    (-(n0 / 2 - 0.5)) * voxel_size(1), ...
    (-(n1 / 2 - 0.5)) * voxel_size(2), ...
    (-(n2 / 2 - 0.5)) * voxel_size(3)]);
clear n0 n1 n2

%% 3. Scanner Geometry / Norm Parameters

% [cry_num_tr  cry_num_ax  mod_num_tr  rsector_num  ring_num]
PET_geom(1,:) = [8 10 6 22 12];

ring_num = 120;
crystal_per_ring = 1056;

min_sector_difference = 2;
crystal_num_per_sector_trans = 48;
min_crystal_difference = crystal_num_per_sector_trans * min_sector_difference + 1; % 97, i.e. > 96
max_ring_difference = ring_num - 1;

num_parts = 22;

%% 4. Parallelproj CUDA Library Configuration

libname = 'libparallelproj';
if ~libisloaded(libname)
    loadlibrary(libname, 'm_parallelproj.h');
end

%% 5. Paths and Files

lm_path = './lmf/';
lut_path = './lut/';
img_path = './img/';
norm_path = './norm/';
initial_mask_path = './initial_and_mask/';

if ~exist(img_path, 'dir')
    mkdir(img_path);
end

% Listmode files
file_lm = [lm_path 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_withI33_20251017_takeAll_511_prompt_halfTime1200s.Cdf'];
file_td = [lm_path 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_withI33_20251017_takeAll_511_prompt_halfTime1200s.timeDiff'];

% LUT
crystal_id_lut_name = 'TBP_noDOI_1axialModule_6p3p4_241028.glut';
PET_lut = [lut_path crystal_id_lut_name];

% Norm factors
NF_ax_file = [norm_path 'NF_ax_TBP_noDOI_1axialModule_1125.mat'];
NF_tr_file = [norm_path 'NF_tr_TBP_noDOI_1axialModule_1125.mat'];

% Initial image / mask
atten_init_file = [initial_mask_path 'Recon_Result_atten_fusion_511_safetyBase_nonNegative_noQP_alpha1_chunks_0p5dose_2mmVoxel_350p138_260426_it20_8mmPostGaussian_masked.img'];
mask_file = [initial_mask_path 'mask_Recon_Result_atten_fusion_511_safetyBase_nonNegative_noQP_alpha1_chunks_0p5dose_2mmVoxel_350p138_260426_it30_6mmPostGaussian_mask0p003.img'];

% Output names
suff_img = '.img';
recon_image_prefix = 'Recon_Result';

sens_initial_file = [img_path 'sens_TB1axial_BrainPhantom_350p138FOV_2mmVoxel_mltrInitial.img'];
activity_initial_file = [img_path 'acti_TB1axial_BrainPhantom_350p138FOV_2mmVoxel_mltrInitial.img'];

%% 6. TOF Parameters: 7 ps bin width

speed_of_light = single(0.3);        % mm / ps

% Coincidence timing resolution FWHM.
time_resolution_ps = single(250);    % ps

% Desired discrete TOF bin width in time.
tofbin_width_ps = single(7);         % ps

% parallelproj expects spatial TOF bin width in mm.
tofbin_width = single(speed_of_light * tofbin_width_ps / 2);  % 1.05 mm

% Gaussian TOF kernel sigma in spatial units.
fwhm_tof_resolution_mm = single(speed_of_light * time_resolution_ps / 2);
sigma_tof = single(fwhm_tof_resolution_mm / (2 * sqrt(2 * log(2))));

tofcenter_offset = single(0);
n_sigmas_TOF = single(3);

lor_dependent_sigma_tof = uint8(0);
lor_dependent_tofcenter_offset = uint8(0);

% Auto-calculate num_tof_bins from timeDiff file.
disp('Scanning timeDiff range to determine num_tof_bins...');

fid_td_scan = fopen(file_td, 'rb');
if fid_td_scan == -1
    error(['Cannot open TD File: ' file_td]);
end

td_scan_chunk = 20000000;
td_min_ps = single(inf);
td_max_ps = single(-inf);

while true
    td_vals = fread(fid_td_scan, td_scan_chunk, 'float=>single');
    if isempty(td_vals)
        break;
    end
    td_min_ps = min(td_min_ps, min(td_vals));
    td_max_ps = max(td_max_ps, max(td_vals));
end

fclose(fid_td_scan);

max_abs_td_ps = max(abs(td_min_ps), abs(td_max_ps));

% Add margin so extreme events do not sit exactly on the boundary bin.
tof_margin_bins = 2;
tof_half_bins = ceil(double(max_abs_td_ps) / double(tofbin_width_ps)) + tof_margin_bins;
num_tof_bins_double = 2 * tof_half_bins + 1;

if num_tof_bins_double > double(intmax('int16'))
    error('num_tof_bins exceeds int16 range required by the CUDA wrapper.');
end

num_tof_bins = int16(num_tof_bins_double);
tof_bin_center_offset = int16(floor(double(num_tof_bins) / 2));

fprintf('\n[TOF parameter summary]\n');
fprintf('  timeDiff min/max        = %.3f / %.3f ps\n', td_min_ps, td_max_ps);
fprintf('  desired TOF bin width   = %.3f ps\n', tofbin_width_ps);
fprintf('  parallelproj width      = %.6f mm\n', tofbin_width);
fprintf('  equivalent width check  = %.6f ps\n', 2 * tofbin_width / speed_of_light);
fprintf('  TOF sigma               = %.6f mm\n', sigma_tof);
fprintf('  TOF FWHM                = %.3f ps, %.6f mm\n', time_resolution_ps, fwhm_tof_resolution_mm);
fprintf('  num_tof_bins            = %d\n', num_tof_bins);
fprintf('  center bin              = %d\n', tof_bin_center_offset);
fprintf('  half time coverage      = %.3f ps\n\n', double(tof_bin_center_offset) * double(tofbin_width_ps));

%% 7. Load LUT and Norm Factors

disp('Loading Crystal LUT...');
crystal_lut = single(readlut(PET_lut));
crystal_lut = crystal_lut(:, 1:3);

disp('Loading Norm Factors...');
NF_ax_data = load(NF_ax_file);
NF_tr_data = load(NF_tr_file);

%% 8. Initial Activity and Attenuation Map

atten_map = single(read_binary_img(atten_init_file));
atten_map = atten_map(:);

mask_map = single(read_binary_img(mask_file));
mask_map = mask_map(:);

% Keep attenuation map inside object support from the start.
atten_map = atten_map .* mask_map;

if initial_load == 0
    image_init_activity = ones(prod(double(img_dim)), 1, 'single');

    sens_all = sensitivity_image_calculate_norm_data_cuda( ...
        libname, crystal_lut, img_dim, voxel_size, img_dim, voxel_size, ...
        atten_map, num_parts, ring_num, crystal_per_ring, ...
        min_crystal_difference, max_ring_difference, ...
        PET_geom, NF_ax_data, NF_tr_data);

    fid_sens = fopen(sens_initial_file, 'wb');
    fwrite(fid_sens, sens_all, 'float32');
    fclose(fid_sens);

    % Initial Activity MLEM
    activity_map = MLEM_parallelproj_TOF_bin_chunk_wAD_cuda( ...
        libname, 10, image_init_activity, ...
        file_lm, file_td, ...
        sens_all, atten_map, ...
        PET_geom, NF_ax_data, NF_tr_data, ...
        crystal_lut, img_origin, voxel_size, img_dim, chunk_size, ...
        tofbin_width, sigma_tof, tofcenter_offset, n_sigmas_TOF, num_tof_bins, ...
        lor_dependent_sigma_tof, lor_dependent_tofcenter_offset);

    fid_activity = fopen(activity_initial_file, 'wb');
    fwrite(fid_activity, activity_map, 'float32');
    fclose(fid_activity);

    clear image_init_activity
else
    sens_all = single(read_binary_img(sens_initial_file));
    sens_all = sens_all(:);

    activity_map = single(read_binary_img(activity_initial_file));
    activity_map = activity_map(:);
end

%% 9. MLAA-TX Iteration Loop

for it = 1:maxit

    disp(['Iteration ' num2str(it) '...']);

    % --- Step 1: Update Attenuation Map ---
    for it_sub = 1:sub_it

        disp(['  Attenuation sub-iteration ' num2str(it_sub) ' / ' num2str(sub_it)]);

        atten_map = ML_attn_sps_parallelproj_lm_norm_chunk_cuda( ...
            libname, atten_map, activity_map, ...
            file_lm, ...
            crystal_lut, PET_geom, NF_ax_data, NF_tr_data, ...
            img_origin, voxel_size, img_dim, chunk_size, ...
            num_parts, ring_num, crystal_per_ring, ...
            min_crystal_difference, max_ring_difference);

        atten_map = atten_map .* mask_map;
    end

    attenuation_out_file = [img_path recon_image_prefix ...
        '_atten_mlaa_TB1axial_BrainPhantom_1200s_350p138FOV_2mmVoxel_IntialMLTR_260707_it' ...
        num2str(it) suff_img];

    fid_atten = fopen(attenuation_out_file, 'wb');
    fwrite(fid_atten, atten_map, 'float32');
    fclose(fid_atten);

    disp('atten map written');

    % --- Step 2: Update Activity Map ---
    sens_all = sensitivity_image_calculate_norm_data_cuda( ...
        libname, crystal_lut, img_dim, voxel_size, img_dim, voxel_size, ...
        atten_map, num_parts, ring_num, crystal_per_ring, ...
        min_crystal_difference, max_ring_difference, ...
        PET_geom, NF_ax_data, NF_tr_data);

    activity_map = MLEM_parallelproj_TOF_bin_chunk_wAD_cuda( ...
        libname, it_local, activity_map, ...
        file_lm, file_td, ...
        sens_all, atten_map, ...
        PET_geom, NF_ax_data, NF_tr_data, ...
        crystal_lut, img_origin, voxel_size, img_dim, chunk_size, ...
        tofbin_width, sigma_tof, tofcenter_offset, n_sigmas_TOF, num_tof_bins, ...
        lor_dependent_sigma_tof, lor_dependent_tofcenter_offset);

    activity_out_file = [img_path recon_image_prefix ...
        '_activ_mlaa_TB1axial_BrainPhantom_1200s_350p138FOV_2mmVoxel_IntialMLTR_260707_it' ...
        num2str(it) suff_img];

    fid_act = fopen(activity_out_file, 'wb');
    fwrite(fid_act, activity_map, 'float32');
    fclose(fid_act);

    disp('activ map written');
    disp(' ');
end
