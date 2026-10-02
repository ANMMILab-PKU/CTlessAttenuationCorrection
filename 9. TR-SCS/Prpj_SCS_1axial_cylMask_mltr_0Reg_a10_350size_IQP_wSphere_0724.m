clc;
clear;

% GPU for non-TOF projection; CPU for TOF projection; all 0 initialization;
initial_load = 0;
disp(mfilename);

%% recon parameters
maxit = 50;
sub_it = 4;
it_local = 1;
chunk_size = 1000000; % 1M events per chunk

% Image Geometry
img_dim = int32([350, 350, 138]);
voxel_size = single([2, 2, 2]);

% Scanner Geometry (For Norm)
PET_geom(1,:) = [8 10 6 22 12];

%% parallel proj
%%%%% Linux cuda %%%%% cuda
libname2 = 'libparallelproj';
if ~libisloaded(libname2), loadlibrary(libname2, 'm_parallelproj.h'); end

%% Configuration
data_path_emits = '/share/home/xzhao/workDir/MLAA_202512/SCS_20260414_GPU/data/data_2607/';
data_path_trans = '/share/home/xzhao/workDir/MLAA_202512/SCS_20260414_GPU/data/data_2607/';
lut_path  = './data/';
img_path  = './img/';
norm_path = './data/';

suff_img = '.img';

% lismode Files
file_em_lm = [data_path_emits 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_NEMA2018_IQP_w50mmBed_wSphere_1800s_takeAll_0716_511_prompt_onlyTrue.Cdf'];
file_em_td = [data_path_emits 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_NEMA2018_IQP_w50mmBed_wSphere_1800s_takeAll_0716_511_prompt_onlyTrue.timeDiff'];

% Unified Histogram Files
file_em     = [data_path_emits 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_NEMA2018_IQP_w50mmBed_wSphere_1800s_takeAll_0716_511_prompt_onlyTrue.histo'];
file_bs_202 = [data_path_trans 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_blankScan_5000s_20251016_takeAll_202_prompt_IQP_1800s_wSphere_2607_unified.histo'];
file_bs_307 = [data_path_trans 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_blankScan_5000s_20251016_takeAll_307_prompt_IQP_1800s_wSphere_2607_unified.histo'];
file_tr_202 = [data_path_trans 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_NEMA2018_IQP_w50mmBed_wSphere_1800s_takeAll_0716_202_prompt_promptAll_unified.histo'];
file_tr_307 = [data_path_trans 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_NEMA2018_IQP_w50mmBed_wSphere_1800s_takeAll_0716_307_prompt_promptAll_unified.histo'];

% 新增additive histo
file_tr_202_add = [data_path_trans 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_NEMA2018_IQP_w50mmBed_wSphere_1800s_takeAll_0716_202_prompt_onlyAdditive_unified.histo'];
file_tr_307_add = [data_path_trans 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_NEMA2018_IQP_w50mmBed_wSphere_1800s_takeAll_0716_307_prompt_onlyAdditive_unified.histo'];

beta_reg = 0;
tx_mode = 3;          % 1=Only202, 2=Only307, 3=Both

alpha_EM_511 = 1;
alpha_TX_202 = 10;
alpha_TX_307 = 10;

%%%%%%%%%%%%%% 1 layer %%%%%%%%%%%%%%%%%%%
ring_num = 120;
crystal_per_ring = 1056;
min_sector_difference = 2;
crystal_num_per_rector_trans = 48;
min_crystal_difference = crystal_num_per_rector_trans*min_sector_difference+1;
max_ring_difference = ring_num - 1;
num_parts = 22;
%%%%%%%%%%%%%% 1 layer %%%%%%%%%%%%%%%%%%%

%% TOF Parameters: 7 ps bin width (matched to MLEM_AC_norm_TOFbin7ps_autoBins_260619)

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

%% Load LUT and Norm
crystal_id_lut_path_name = 'TBP_noDOI_1axialModule_6p3p4_241028.glut';
PET_lut = [lut_path crystal_id_lut_path_name];
crystal_lut = readlut(PET_lut);
crystal_lut = single(crystal_lut(:,1:3));

recon_image_path_name_all = 'Recon_Result';
NF_ax_file = [norm_path 'NF_ax_TBP_noDOI_1axialModule_1125.mat'];
NF_tr_file = [norm_path 'NF_tr_TBP_noDOI_1axialModule_1125.mat'];

disp('Loading Norm Factors into memory...');
NF_ax_data = load(NF_ax_file);
NF_tr_data = load(NF_tr_file);

%% Image origin setup
n0 = single(img_dim(1)); n1 = single(img_dim(2)); n2 = single(img_dim(3));
img_origin_activity = single([(-(n0 / 2 - 0.5)) * voxel_size(1), (-(n1 / 2 - 0.5)) * voxel_size(2), (-(n2 / 2 - 0.5)) * voxel_size(3)]);
clear n0 n1 n2

%% 3. Initialization
% 加载初始或掩码 atten_map
atten_map = read_binary_img('/share/home/xzhao/workDir/MLAA_202512/SCS_20260414_GPU/data/Recon_Result_atten_fusion_511_200QP_alpha1_additCorr_chunks_2mmVoxel_350p138_260723_it17_6mmPostGaussian_masked.img');
atten_map = atten_map(:);

mask_map = read_binary_img('/share/home/xzhao/workDir/MLAA_202512/SCS_20260414_GPU/data/mask_cylinder350_350p138_260521.img');
mask_map = mask_map(:);


if initial_load == 0
    image_init_activity = single(ones(img_dim));
    image_init_activity = image_init_activity(:);
    sens_all = sensitivity_image_calculate_norm_data_cuda(libname2, crystal_lut, img_dim, voxel_size, img_dim, voxel_size, atten_map, num_parts, ring_num, crystal_per_ring, min_crystal_difference, max_ring_difference, PET_geom, NF_ax_data, NF_tr_data);

    fid_sens = fopen('sens_TB1axial_IQPw50mmBed_wSphere_260717_350p138FOV_2mmVoxel_mltrInitial.img', 'wb');
    fwrite(fid_sens, sens_all, 'float32'); fclose(fid_sens);

    activity_map = MLEM_parallelproj_TOF_bin_chunk_wAD_cuda(...
        libname2, 15, image_init_activity, file_em_lm, file_em_td, sens_all, ...
        atten_map, PET_geom, NF_ax_data, NF_tr_data, crystal_lut, ...
        img_origin_activity, voxel_size, img_dim, chunk_size, ...
        tofbin_width, sigma_tof, tofcenter_offset, n_sigmas_TOF, num_tof_bins, ...
        lor_dependent_sigma_tof, lor_dependent_tofcenter_offset);

    fid_activity = fopen('acti_TB1axial_IQPw50mmBed_wSphere_260717_350p138FOV_2mmVoxel_mltrInitial.img', 'wb');
    fwrite(fid_activity, activity_map, 'float32'); fclose(fid_activity);
    clear image_init_activity
elseif initial_load == 1
    sens_all     = read_binary_img('sens_TB1axial_IQPw50mmBed_wSphere_260717_350p138FOV_2mmVoxel_mltrInitial.img'); sens_all = sens_all(:);
    activity_map = read_binary_img('acti_TB1axial_IQPw50mmBed_wSphere_260717_350p138FOV_2mmVoxel_mltrInitial.img'); activity_map = activity_map(:);
end

% =========================================================================
% === [新增] Pass 0: Initialization (SCS Algorithm)                     ===
% =========================================================================
fprintf('\n================================================================\n');
fprintf('=== Pass 0: Initialization (SCS Algorithm)                   ===\n');
fprintf('================================================================\n\n');

params_SCS = struct();
% 【已移除硬编码的eta_202和eta_307】
params_SCS.alpha_EM  = alpha_EM_511;
params_SCS.alpha_202 = alpha_TX_202;
params_SCS.alpha_307 = alpha_TX_307;
params_SCS.beta_reg  = beta_reg;
params_SCS.tx_mode   = tx_mode;
params_SCS.scale_factor = (5000 * 7) / (1800); % 转换 BS 到 TX 尺度的因子
params_SCS.chunk_size = chunk_size;
params_SCS.file_tr_202_add = file_tr_202_add;
params_SCS.file_tr_307_add = file_tr_307_add;

psi_geom = struct();
psi_geom.libname = libname2;
psi_geom.crystal_lut = crystal_lut;
psi_geom.img_dim = img_dim;
psi_geom.voxel_size = voxel_size;
psi_geom.img_origin = img_origin_activity;
psi_geom.PET_geom = PET_geom;
psi_geom.NF_ax_data = NF_ax_data;
psi_geom.NF_tr_data = NF_tr_data;
psi_geom.num_parts = num_parts;
psi_geom.ring_num = ring_num;
psi_geom.crystal_per_ring = crystal_per_ring;
psi_geom.min_crystal_difference = min_crystal_difference;
psi_geom.max_ring_difference = max_ring_difference;

scs_state = struct();
scs_state.lambda_cubic = 1.0;
scs_state.M_j_analytical = []; 
scs_state.M_j_previous = [];

diagnostics = struct('grad_norm', zeros(maxit, sub_it), 'step_norm', zeros(maxit, sub_it), 'hess_pos_ratio', zeros(maxit, sub_it), 'lambda', zeros(maxit, sub_it));
fprintf('Pass 0 Complete.\n\n');


[FI, alpha_rec] = calc_fisher_information_SCS(atten_map, activity_map, file_em, file_tr_202, file_bs_202, file_tr_307, file_bs_307, params_SCS, psi_geom);

% 应用推荐权重（可乘以 gamma 因子微调 TX 相对权重）
% gamma = 2.0;  % gamma=1 为等 FI 贡献，>1 偏向 TX
% params_SCS.alpha_EM  = alpha_rec.alpha_EM;
% params_SCS.alpha_202 = alpha_rec.alpha_202 * gamma;
% params_SCS.alpha_307 = alpha_rec.alpha_307 * gamma;

% =========================================================================
% 4. MLAA-TX Iteration Loop (Alternating Optimization)
% =========================================================================
for it = 1:maxit
    disp(['Iteration ' num2str(it) '...']);

    % --- Step 1: Update Attenuation (SCS) ---
    tic;
    for it_sub = 1:sub_it
        fprintf('  Sub-iteration %d/%d (SCS Attenuation Update)...\n', it_sub, sub_it);

        [atten_map, scs_state, dbg_info] = MLAA_TX_SCS_update_MEP_260518(atten_map, activity_map, file_em, file_tr_202, file_bs_202, file_tr_307, file_bs_307, params_SCS, psi_geom, scs_state, it, it_sub);
     	atten_map = atten_map .* mask_map;

        diagnostics.grad_norm(it, it_sub) = dbg_info.grad_norm;
        diagnostics.step_norm(it, it_sub) = dbg_info.step_norm;
        diagnostics.hess_pos_ratio(it, it_sub) = dbg_info.hess_pos_ratio;
        diagnostics.lambda(it, it_sub) = scs_state.lambda_cubic;

        fprintf('    ||g||=%.4e, HessPos=%.2f%%, SaddleEscapes=%d, lambda=%.4f\n', dbg_info.grad_norm, 100*dbg_info.hess_pos_ratio, dbg_info.saddle_escapes, scs_state.lambda_cubic);

    end
    t_mu = toc;
    disp(['  Mu-map update done in ' num2str(t_mu) 's']);

    % Save Mu-map
    fid_write_attenuation = fopen([img_path recon_image_path_name_all '_atten_w50mmBed_wSphere_SCS_TBPET1axial_mltrInital_350p138FOV_cylMask_0Reg_alpha10_1800s_additCorr_3MEP_0724_it' num2str(it) suff_img], 'wb');

    fwrite(fid_write_attenuation, atten_map, 'float32'); fclose(fid_write_attenuation);
    disp('atten map writen');

    % --- Step 2: Update Activity (Standard MLEM) ---
    tic;
    sens_all = sensitivity_image_calculate_norm_data_cuda(libname2, crystal_lut, img_dim, voxel_size, img_dim, voxel_size,  atten_map, num_parts, ring_num, crystal_per_ring, min_crystal_difference, max_ring_difference, PET_geom, NF_ax_data, NF_tr_data);

    activity_map = MLEM_parallelproj_TOF_bin_chunk_wAD_cuda(...
        libname2, it_local, activity_map, file_em_lm, file_em_td, sens_all, ...
        atten_map, PET_geom, NF_ax_data, NF_tr_data, crystal_lut, ...
        img_origin_activity, voxel_size, img_dim, chunk_size, ...
        tofbin_width, sigma_tof, tofcenter_offset, n_sigmas_TOF, num_tof_bins, ...
        lor_dependent_sigma_tof, lor_dependent_tofcenter_offset);

    activity_map = activity_map .* mask_map;

    t_ac = toc;
    disp(['  activity-map update done in ' num2str(t_ac) 's']);

    fid_write_activity = fopen([img_path recon_image_path_name_all '_activ_w50mmBed_wSphere_SCS_TBPET1axial_mltrInital_350p138FOV_cylMask_0Reg_alpha10_1800s_additCorr_3MEP_0724_it' num2str(it) suff_img], 'wb');
    fwrite(fid_write_activity, activity_map, 'float32'); fclose(fid_write_activity);
    disp('activ map writen');
    disp('    ');
end
