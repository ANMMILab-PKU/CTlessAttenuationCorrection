function Ratio_factor_image = proj_back_parallelproj_cuda(libname, xstart, xend, image_recon, img_origin, voxel_size, Ratio_factor_sino, nlors, img_dim)
%--------------------------------------------------------------------------
% back projection 
%--------------------------------------------------------------------------

% 准备数据...
image_size_xy = img_dim(1);
image_size_z = img_dim(3);

image_recon = zeros(size(image_recon));
image_recon = single(image_recon(:));

Ratio_factor_sino = single(Ratio_factor_sino(:));

% 新增 CUDA 参数
device_id = int32(0);
threads = int32(64);

nlors     = uint64(nlors);

% 调用 CUDA 库，末尾增加 device_id 和 threads
[~, ~, Ratio_factor_image, ~, ~, ~, ~] = calllib(libname, 'joseph3d_back', xstart, xend, image_recon, img_origin, voxel_size, Ratio_factor_sino, nlors, img_dim, device_id, threads);

Ratio_factor_image = reshape(Ratio_factor_image, [image_size_z, image_size_xy, image_size_xy]);
Ratio_factor_image = permute(Ratio_factor_image, [3, 2, 1]);

Ratio_factor_image = double(Ratio_factor_image(:));
