% clc
clear

initial_load = 0;

%% recon parameters
maxit = 100;
chunk_size = 500000;

%% image parameters
img_dim = int32([350, 350, 138]);
voxel_size = single([2, 2, 2]);

%%%%%%%%%%%%%% 1 layer %%%%%%%%%%%%%%%%%%%
ring_num = 120;
crystal_per_ring = 1056;
min_sector_difference = 2;
crystal_num_per_rector_trans = 48;
min_crystal_diffenrence = crystal_num_per_rector_trans * min_sector_difference + 1;
max_ring_difference = ring_num - 1;
num_parts = 22;
diemeter_scanner = 410.9 * 2; % 完美修正：这里现在是真实的物理直径！

% Paper setting is alpha = 1.5.
% If the update is too aggressive for your data scale, try alpha_step = 1.0.
alpha_step = 1.5;
%%%%%%%%%%%%%% 1 layer %%%%%%%%%%%%%%%%%%%

%% parallel proj
%%%%% Linux %%%%%
libname = 'libparallelproj_c';
tf = libisloaded(libname);
if ~tf
    loadlibrary(libname, 'parallelproj_c.h');
end
%%%%% Linux %%%%%

%% data paths
data_path_emits = '/share/home/xzhao/workDir/MLAA_202512/MLTR_20260419_GPU/data/';
data_path_trans = '/share/home/xzhao/workDir/MLAA_202512/MLTR_20260419_GPU/data/';
lut_path  = './data/';
img_path  = './img/';

if ~exist(img_path, 'dir')
    mkdir(img_path);
end

% Unified Histogram Files
file_bs_202 = [data_path_trans 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_blankScan_5000s_20251016_takeAll_202_prompt_IQP_unified.histo'];
file_bs_307 = [data_path_trans 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_blankScan_5000s_20251016_takeAll_307_prompt_IQP_unified.histo'];

file_tr_202 = [data_path_trans 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_NEMA2018_IQP_w50mmBed_1800s_takeAll_0619_202_prompt_onlyAll_unified.histo'];
file_tr_307 = [data_path_trans 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_NEMA2018_IQP_w50mmBed_1800s_takeAll_0619_307_prompt_onlyAll_unified.histo'];

% Additive histograms
file_tr_202_add = [data_path_trans 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_NEMA2018_IQP_w50mmBed_1800s_takeAll_0619_202_prompt_onlyNotTrue_unified.histo'];
file_tr_307_add = [data_path_trans 'coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_NEMA2018_IQP_w50mmBed_1800s_takeAll_0619_307_prompt_onlyNotTrue_unified.histo'];

crystal_id_lut_path_name = 'TBP_noDOI_1axialModule_6p3p4_241028.glut';
suff_img = '.img';
recon_image_path_name_all = 'Recon_Result';

crystal_lut = readlut([lut_path crystal_id_lut_path_name]);
crystal_lut = single(crystal_lut(:, 1:3));

%% Time scale factor
scale_factor = single((5000 * 7) / 1800);

%% QP regularization
beta_qp = single(0);
factor = single(1.0);
use_factor_ramp = false;
factor_ramp_iters = 10;

beta_tag = strrep(num2str(double(beta_qp)), '.', 'p');
factor_tag = strrep(num2str(double(factor)), '.', 'p');


%% =========================================================
%% 新增：生成 Cylinder Mask (在循环外只计算一次，极大地提升效率)
%% =========================================================
disp('Generating Cylinder FOV Mask...');
% Mask 直径 = 正方形边长
fov_diameter = single(img_dim(1)) * voxel_size(1); % 300 * 2 = 600 mm
fov_radius = fov_diameter / single(2.0); % 半径 300 mm

% 生成二维中心对称坐标网格
x_coords = (-(single(img_dim(1))/2 - 0.5)) * voxel_size(1) : voxel_size(1) : (single(img_dim(1))/2 - 0.5) * voxel_size(1);
y_coords = (-(single(img_dim(2))/2 - 0.5)) * voxel_size(2) : voxel_size(2) : (single(img_dim(2))/2 - 0.5) * voxel_size(2);
[X_grid, Y_grid] = ndgrid(x_coords, y_coords);

% 计算半径的平方并生成 Mask
R_sq = single(X_grid.^2 + Y_grid.^2);
mask_2d = R_sq <= (fov_radius^2);

% 将 2D Mask 扩展到 3D
mask_3d = repmat(mask_2d, [1, 1, img_dim(3)]);

% 转换为一维向量 (由于在调用函数前将图像转为一维，这里配合拉直)
mask_1d = mask_3d(:); 
disp('Mask generated successfully.');


%% initialization
if initial_load == 0
    image_init_attenuation = single(ones(img_dim));
    image_init_attenuation = image_init_attenuation .* single(0.0096); % water at 511 keV, 1/mm
    
    % 在初始化阶段，就把 FOV 外部的水设置为空，防止前向投影产生不必要的衰减
    image_init_attenuation(~mask_3d) = 0; 
    
    image_recon_attenuation = image_init_attenuation(:);
    clear image_init_attenuation
end

image_recon_attenuation_202 = image_recon_attenuation;
image_recon_attenuation_307 = image_recon_attenuation;
clear image_recon_attenuation

%% iterative loop
t1 = tic;

for it = 1:maxit
    disp(['iteration: ' num2str(it)]);

    if use_factor_ramp
        factor_it = factor .* min(single(1), single(it) ./ single(factor_ramp_iters));
    else
        factor_it = factor;
    end

    disp(['  Current QP factor_it = ' num2str(double(factor_it)) ...
          ', beta_eff = ' num2str(double(beta_qp .* factor_it))]);

    % -----------------------------------------------------------
    % 202 Reconstruction with additive correction + QP
    % -----------------------------------------------------------
    disp('  Processing 202 dataset...');
    image_recon_attenuation_202 = MLTR_parallelproj_histogram_chunk_wAD_QP_260623_denomfix( ...
        libname, image_recon_attenuation_202, file_bs_202, file_tr_202, file_tr_202_add, ...
        crystal_lut, img_dim, voxel_size, diemeter_scanner, alpha_step, ...
        chunk_size, scale_factor, beta_qp, factor_it);

    % 新增：更新完成后施加 Mask 截断
    image_recon_attenuation_202(~mask_1d) = 0;

    fid_write_202 = fopen([img_path recon_image_path_name_all '_atten_202_additCorr_QP_beta' beta_tag '_factor' factor_tag '_260624_350p138_it' num2str(it) suff_img], 'wb');
    fwrite(fid_write_202, image_recon_attenuation_202, 'float32');
    fclose(fid_write_202);

    % -----------------------------------------------------------
    % 307 Reconstruction with additive correction + QP
    % -----------------------------------------------------------
    disp('  Processing 307 dataset...');
    image_recon_attenuation_307 = MLTR_parallelproj_histogram_chunk_wAD_QP_260623_denomfix( ...
        libname, image_recon_attenuation_307, file_bs_307, file_tr_307, file_tr_307_add, crystal_lut, img_dim, voxel_size, diemeter_scanner, alpha_step, ...
        chunk_size, scale_factor, beta_qp, factor_it);

    % 新增：更新完成后施加 Mask 截断
    image_recon_attenuation_307(~mask_1d) = 0;

    fid_write_307 = fopen([img_path recon_image_path_name_all '_atten_307_additCorr_QP_beta' beta_tag '_factor' factor_tag '_260624_350p138_it' num2str(it) suff_img], 'wb');
    fwrite(fid_write_307, image_recon_attenuation_307, 'float32');
    fclose(fid_write_307);

    disp('  Atten maps written for this iteration.');
    disp('    ');
end

t2 = toc(t1);
disp(['totally time consumption: ' num2str(t2) ' s']);
