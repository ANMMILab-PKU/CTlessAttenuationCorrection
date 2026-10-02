function sino_fwd = proj_forw_parallelproj_cuda(libname, xstart, xend, image_recon, img_origin, voxel_size, sino_init, nlors, img_dim)
%--------------------------------------------------------------------------
% forward projection
%--------------------------------------------------------------------------

% 准备数据...
image_size_xy = img_dim(1);
image_size_z = img_dim(3);

image_recon = reshape(image_recon, [image_size_xy, image_size_xy, image_size_z]);
image_recon = permute(image_recon, [3, 2, 1]);

image_recon = single(image_recon(:));

sino_init=zeros(size(sino_init));
sino_init = single(sino_init(:));

% 新增 CUDA 参数
device_id = int32(0);
threads = int32(64);

nlors     = uint64(nlors);   % size_t

% 调用 CUDA 库，末尾增加 device_id 和 threads
[~, ~, ~, ~, ~, sino_fwd, ~] = calllib(libname, 'joseph3d_fwd', xstart, xend, image_recon, img_origin, voxel_size, sino_init, nlors, img_dim, device_id, threads);

sino_fwd = double(sino_fwd(:));