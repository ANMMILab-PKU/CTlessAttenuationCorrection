function sens_all = sensitivity_image_calculate_norm_data_cuda(libname, crystal_lut, ...
    img_dim_sens, voxel_size_sens, img_dim_attenuation, voxel_size_miumap, ...
    atten_map, num_parts, ring_num, crystal_per_ring, min_crystal_diffenrence, max_ring_difference, ...
    PET_geom, NF_ax_data, NF_tr_data) % [新增参数] 接收几何和归一化数据

% SENSITIVITY_IMAGE_CALCULATE_NORM
% 修改版：不再在内部 load 文件，而是接收预加载的数据结构。

% [删除] 以前的硬编码路径和 load 逻辑
% PET_lut = '/path/to/glut'; 
% PET_geom(1,:) = [8 10 6 22 12];  <-- 几何参数现在从外部传入，保证一致性
% NF_ax = '/path/to/NF_ax.mat'; 
% NF_tr = '/path/to/NF_tr.mat';

%% 参数初始化
use_psf = false;
FWHM = [1, 1, 1];
psf_type = 'version2';
nbSigmas = 3;

% 预分配内存 
sens_all = single(zeros(prod([img_dim_sens(1), img_dim_sens(2), img_dim_sens(3)]), 1));       

for i_part = 1 : num_parts

    % 计算有效 LOR
    sino_id = Get_valid_lor_calculate_sens_single_layer(ring_num, crystal_per_ring, min_crystal_diffenrence, max_ring_difference, num_parts, i_part);
    
    % 利用传入的 crystal_lut 数组获取坐标
    coor_sinoid1 = crystal_lut(sino_id(:, 1), :);
    coor_sinoid2 = crystal_lut(sino_id(:, 2), :);

    xstart = coor_sinoid1';
    xend = coor_sinoid2';
    xstart = single(xstart(:));
    xend = single(xend(:));
    clear coor_sinoid1 coor_sinoid2

    % 维度转换
    img_dim_sens    = int32(img_dim_sens);
    img_dim_miumap  = int32(img_dim_attenuation);

    n0_miumap = single(img_dim_miumap(1));
    n1_miumap = single(img_dim_miumap(2));
    n2_miumap = single(img_dim_miumap(3));

    img_origin_sens     = single([(-(single(img_dim_sens(1)) / 2 - 0.5)) * voxel_size_sens(1), (-(single(img_dim_sens(2)) / 2 - 0.5)) * voxel_size_sens(2), (-(single(img_dim_sens(3)) / 2 - 0.5)) * voxel_size_sens(3)]);
    img_origin_miumap   = single([(-(n0_miumap / 2 - 0.5)) * voxel_size_miumap(1), (-(n1_miumap / 2 - 0.5)) * voxel_size_miumap(2), (-(n2_miumap / 2 - 0.5)) * voxel_size_miumap(3)]);

    nlors = size(xstart, 1) / 3;

    sino_init = single(ones((size(xstart, 1)) / 3, 1));
    image_init_sens = single(ones(img_dim_sens));
    image_init_sens = image_init_sens(:);

    % 衰减因子计算 (Forward Proj)
    sino_fwd = proj_forw_parallelproj_cuda(libname, xstart, xend, atten_map, img_origin_miumap, voxel_size_miumap, sino_init, nlors, img_dim_miumap);
    
    % line_integral = sino_fwd ./ 10; % cm-1 * mm ---> cm-1 * cm
    line_integral = sino_fwd; % mm-1
    attenuation_factor = exp((-1) * line_integral);    
    sino_init = sino_init .* attenuation_factor;
    
    % [修改] 调用 NF_generate
    LORs_frac = sino_id - 1;
    
    % 注意：这里直接传入 crystal_lut (数组) 和 NF 数据结构
    normalization_factor = NF_generate_data(LORs_frac, crystal_lut, PET_geom, NF_ax_data, NF_tr_data);
    
    sino_init = sino_init .* (1./normalization_factor);
    sino_init(isnan(sino_init) | isinf(sino_init)) = 0;
    
    % 反投影 (Back Proj)
    sens = proj_back_parallelproj_cuda(libname, xstart, xend, sens_all, img_origin_sens, voxel_size_sens, sino_init, nlors, img_dim_sens);

    sens_all = sens_all + sens;
    clear sens
end

end
