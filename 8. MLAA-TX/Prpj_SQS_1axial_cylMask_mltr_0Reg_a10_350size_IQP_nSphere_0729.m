clc;
clear;

% GPU for non-TOF projection; CPU for TOF projection; all 0 initialization;

initial_load = 1;

disp(mfilename);

%% 1. Recon parameters

maxit = 50;
resume_from_it = 20; % Load the completed it12 checkpoint and continue at it13
sub_it = 4;

it_local = 1;

chunk_size = 10000000; % 1M events per chunk

% Image Geometry
image_size_x = 350;
image_size_y = 350;
image_size_z = 138;

img_dim = int32([image_size_x, image_size_y, image_size_z]);
clear image_size_x image_size_y image_size_z

% voxel size
voxel_size_x = 2;
voxel_size_y = 2;
voxel_size_z = 2;

voxel_size = single([voxel_size_x, voxel_size_y, voxel_size_z]);
clear voxel_size_x voxel_size_y voxel_size_z

% Scanner Geometry (For Norm)
PET_geom(1,:) = [8 10 6 22 12];

%% 2. Parallel proj lib configuration
%%%%% Linux cuda %%%%%
libname2 = 'libparallelproj';
if ~libisloaded(libname2), loadlibrary(libname2, 'm_parallelproj.h'); end

%% 3. Configuration & Paths

% Paths
data_path_emits = '/share/home/xzhao/workDir/MLAA_202512/SCS_20260414_GPU/data/data_2607/';
data_path_trans = '/share/home/xzhao/workDir/MLAA_202512/SCS_20260414_GPU/data/data_2607/';
lut_path  = './data/';
img_path  = './img/';
norm_path = './data/';

% lismode Files
file_em_lm = [data_path_emits 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_NEMA2018_IQP_w50mmBed_noSphere_1800s_takeAll_0717_511_prompt_onlyTrue.Cdf'];
file_em_td = [data_path_emits 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_NEMA2018_IQP_w50mmBed_noSphere_1800s_takeAll_0717_511_prompt_onlyTrue.timeDiff'];

% Unified Histogram Files
file_em     = [data_path_emits 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_NEMA2018_IQP_w50mmBed_noSphere_1800s_takeAll_0717_511_prompt_onlyTrue.histo'];
file_bs_202 = [data_path_trans 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_blankScan_5000s_20251016_takeAll_202_prompt_IQP_1800s_noSphere_2607_unified.histo'];
file_bs_307 = [data_path_trans 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_blankScan_5000s_20251016_takeAll_307_prompt_IQP_1800s_noSphere_2607_unified.histo'];
file_tr_202 = [data_path_trans 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_NEMA2018_IQP_w50mmBed_noSphere_1800s_takeAll_0717_202_prompt_promptAll_unified.histo'];
file_tr_307 = [data_path_trans 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_NEMA2018_IQP_w50mmBed_noSphere_1800s_takeAll_0717_307_prompt_promptAll_unified.histo'];

% 新增additive histo
file_tr_202_add = [data_path_trans 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_NEMA2018_IQP_w50mmBed_noSphere_1800s_takeAll_0717_202_prompt_onlyAdditive_unified.histo'];
file_tr_307_add = [data_path_trans 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_NEMA2018_IQP_w50mmBed_noSphere_1800s_takeAll_0717_307_prompt_onlyAdditive_unified.histo'];

% Energy Scaling Factors (511keV 归一化)
% eta_202 = 0.135 / 0.096; % ~1.40
% eta_307 = 0.118 / 0.096; % ~1.23

% scale_factor = (10000 * 10) / (240 * 10); % Blank/Phantom time scale  % ≈41.67
% scale_factor =  (5000 * 7) / (600 * 2); % ≈ 29.17
scale_factor =  (5000 * 7) / (1800); %


tx_mode = 3; % 1=Only202, 2=Only307, 3=Both

alpha_EM_511 = 1;
alpha_TX_202 = 10;
alpha_TX_307 = 10;
beta_reg = 0; % Regularization weight

% Scanner Geometry 1-axial
ring_num = 120;
crystal_per_ring = 1056;
min_sector_difference = 2;
crystal_num_per_rector_trans = 48;
min_crystal_difference = crystal_num_per_rector_trans*min_sector_difference+1;
max_ring_difference = ring_num - 1;
num_parts = 22;

%% 4. TOF Parameters: 7 ps bin width (matched to MLEM_AC_norm_TOFbin7ps_autoBins_260619)

speed_of_light = single(0.3);        % mm / ps

% Coincidence timing resolution FWHM (not the TOF bin width).
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

% --- Auto-calculate num_tof_bins from timeDiff file ---
disp('Scanning timeDiff range to determine num_tof_bins...');

fid_td_scan = fopen(file_em_td, 'rb');
if fid_td_scan == -1
    error(['Cannot open TD File: ' file_em_td]);
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

% Add a small margin so that the most extreme events do not sit
% exactly on the boundary bin.
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

%% 5. Load LUT & Norm Factors

crystal_id_lut_path_name = 'TBP_noDOI_1axialModule_6p3p4_241028.glut';
PET_lut = [lut_path crystal_id_lut_path_name];
crystal_lut = single(readlut(PET_lut));
crystal_lut = crystal_lut(:,1:3);

suff_img = '.img';
recon_image_path_name_all = 'Recon_Result';

NF_ax_file = [norm_path 'NF_ax_TBP_noDOI_1axialModule_1125.mat'];
NF_tr_file = [norm_path 'NF_tr_TBP_noDOI_1axialModule_1125.mat'];

disp('Loading Norm Factors...');
NF_ax_data = load(NF_ax_file);
NF_tr_data = load(NF_tr_file);

%% 6. Initial Activity & Attenuation Map

% tofcenter_offset already set in TOF Parameters section above
n0 = single(img_dim(1)); n1 = single(img_dim(2)); n2 = single(img_dim(3));
img_origin_activity = single([(-(n0 / 2 - 0.5)) * voxel_size(1), (-(n1 / 2 - 0.5)) * voxel_size(2), (-(n2 / 2 - 0.5)) * voxel_size(3)]);
clear n0 n1 n2

mask_map = read_binary_img('/share/home/xzhao/workDir/MLAA_202512/SCS_20260414_GPU/data/mask_cylinder350_350p138_260521.img');
mask_map = mask_map(:);

if resume_from_it > 0
    % A complete global checkpoint consists of both the attenuation and
    % activity images written at the end of the same iteration.
    atten_checkpoint = fullfile(img_path, [recon_image_path_name_all '_atten_w50mmBed_noSphere_SQS_TBPET1axial_mltrInital_350p138FOV_cylMask_0Reg_alpha10_1800s_additCorr_3MEP_0724_it'  num2str(resume_from_it) suff_img]);
    activity_checkpoint = fullfile(img_path, [recon_image_path_name_all   '_activ_w50mmBed_noSphere_SQS_TBPET1axial_mltrInital_350p138FOV_cylMask_0Reg_alpha10_1800s_additCorr_3MEP_0724_it'   num2str(resume_from_it) suff_img]);

    if exist(atten_checkpoint, 'file') ~= 2
        error('Cannot find attenuation checkpoint: %s', atten_checkpoint);
    end
    if exist(activity_checkpoint, 'file') ~= 2
        error('Cannot find activity checkpoint: %s', activity_checkpoint);
    end

    expected_voxels = prod(double(img_dim));
    expected_bytes = 4 * expected_voxels; % float32 raw image
    atten_info = dir(atten_checkpoint);
    activity_info = dir(activity_checkpoint);
    if atten_info.bytes ~= expected_bytes
        error('Incomplete attenuation checkpoint: expected %d bytes, found %d.', expected_bytes, atten_info.bytes);
    end
    if activity_info.bytes ~= expected_bytes
        error('Incomplete activity checkpoint: expected %d bytes, found %d.', expected_bytes, activity_info.bytes);
    end

    atten_map = read_binary_img(atten_checkpoint);
    atten_map = single(atten_map(:));
    activity_map = read_binary_img(activity_checkpoint);
    activity_map = single(activity_map(:));

    if numel(atten_map) ~= expected_voxels || numel(activity_map) ~= expected_voxels
        error('Checkpoint dimensions do not match img_dim.');
    end

    fprintf('Loaded completed global iteration %d; continuing at iteration %d.\n', resume_from_it, resume_from_it + 1);
else
    atten_map = read_binary_img('/share/home/xzhao/workDir/MLAA_202512/SCS_20260414_GPU/data/Recon_Result_atten_fusion_511_200QP_alpha1_additCorr_chunks_2mmVoxel_350p138_260724_it17_6mmPostGaussian_masked.img');
    atten_map = atten_map(:);

    if initial_load == 0
        image_init_activity = single(ones(img_dim));
        image_init_activity = image_init_activity(:);
        sens_all = sensitivity_image_calculate_norm_data_cuda(libname2, crystal_lut, img_dim, voxel_size, img_dim, voxel_size, atten_map, num_parts, ring_num, crystal_per_ring, min_crystal_difference, max_ring_difference, PET_geom, NF_ax_data, NF_tr_data);

        fid_sens = fopen('sens_TB1axial_IQPw50mmBed_noSphere_260717_350p138FOV_2mmVoxel_mltrInitial.img', 'wb');
        fwrite(fid_sens, sens_all, 'float32');
        fclose(fid_sens);

        % 初次 Activity MLEM
        activity_map = MLEM_parallelproj_TOF_bin_chunk_wAD_cuda(...
            libname2, 15, image_init_activity, file_em_lm, file_em_td, sens_all, ...
            atten_map, PET_geom, NF_ax_data, NF_tr_data, crystal_lut, ...
            img_origin_activity, voxel_size, img_dim, chunk_size, ...
            tofbin_width, sigma_tof, tofcenter_offset, n_sigmas_TOF, num_tof_bins, ...
            lor_dependent_sigma_tof, lor_dependent_tofcenter_offset);

        fid_activity = fopen('acti_TB1axial_IQPw50mmBed_noSphere_260717_350p138FOV_2mmVoxel_mltrInitial.img', 'wb');
        fwrite(fid_activity, activity_map, 'float32');
        fclose(fid_activity);

        clear image_init_activity
    else
        sens_all = read_binary_img('sens_TB1axial_IQPw50mmBed_noSphere_260717_350p138FOV_2mmVoxel_mltrInitial.img');
        sens_all = sens_all(:);
        activity_map = read_binary_img('acti_TB1axial_IQPw50mmBed_noSphere_260717_350p138FOV_2mmVoxel_mltrInitial.img');
        activity_map = activity_map(:);
    end
end

start_it = resume_from_it + 1;

%% 7. MLAA-TX Iteration Loop

for it = start_it:maxit

    disp(['Iteration ' num2str(it) '...']);

    % --- Step 1: Update Attenuation (MLAA-TX SQS Unified) ---
    tic;
    for it_sub = 1:sub_it

        activity_map = activity_map .* mask_map;
        atten_map = MLAA_TX_parallelproj_histo_chunk_data_wAD_LCcul_3MEP_scale(...
            libname2, atten_map, activity_map, file_em, ...
            file_tr_202, file_bs_202, file_tr_202_add, ...
            file_tr_307, file_bs_307, file_tr_307_add, ...
            PET_lut, crystal_lut, img_dim, voxel_size, ...
            alpha_EM_511, alpha_TX_202, alpha_TX_307, ...
            beta_reg, chunk_size, tx_mode, ...
            PET_geom, NF_ax_data, NF_tr_data, ...
            num_parts, ring_num, crystal_per_ring, ...
            min_crystal_difference, max_ring_difference, scale_factor);
        atten_map = atten_map .* mask_map;

    end
    t_mu = toc;
    disp(['  Mu-map update done in ' num2str(t_mu) 's']);

    % Save Mu-map
    fid_mu = fopen([img_path recon_image_path_name_all '_atten_w50mmBed_noSphere_SQS_TBPET1axial_mltrInital_350p138FOV_cylMask_0Reg_alpha10_1800s_additCorr_3MEP_0724_it' num2str(it) suff_img], 'wb');


    fwrite(fid_mu, atten_map, 'float32');
    fclose(fid_mu);

    % --- Step 2: Update Activity (Standard MLEM) ---
    tic;
    % 根据新 Mu-map 重新计算敏感度
    sens_all = sensitivity_image_calculate_norm_data_cuda(libname2, crystal_lut, img_dim, voxel_size, img_dim, voxel_size, atten_map, num_parts, ring_num, crystal_per_ring, min_crystal_difference, max_ring_difference, PET_geom, NF_ax_data, NF_tr_data);

    activity_map = MLEM_parallelproj_TOF_bin_chunk_wAD_cuda(...
        libname2, it_local, activity_map, file_em_lm, file_em_td, sens_all, ...
        atten_map, PET_geom, NF_ax_data, NF_tr_data, crystal_lut, ...
        img_origin_activity, voxel_size, img_dim, chunk_size, ...
        tofbin_width, sigma_tof, tofcenter_offset, n_sigmas_TOF, num_tof_bins, ...
        lor_dependent_sigma_tof, lor_dependent_tofcenter_offset);
    activity_map = activity_map .* mask_map;

    t_ac = toc;
    disp(['  Activity-map update done in ' num2str(t_ac) 's']);

    % Save Activity-map
    fid_act = fopen([img_path recon_image_path_name_all '_activ_w50mmBed_noSphere_SQS_TBPET1axial_mltrInital_350p138FOV_cylMask_0Reg_alpha10_1800s_additCorr_3MEP_0724_it' num2str(it) suff_img], 'wb');

    fwrite(fid_act, activity_map, 'float32');
    fclose(fid_act);

    disp(' ');
end
