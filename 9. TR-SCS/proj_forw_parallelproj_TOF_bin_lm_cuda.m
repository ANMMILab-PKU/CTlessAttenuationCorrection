function sino_fwd = proj_forw_parallelproj_TOF_bin_lm_cuda(libname, xstart, xend, image_recon, img_origin, voxel_size, sino_init, nlors, img_dim, tofbin_width, sigma_tof, tofcenter_offset, n_sigmas_TOF, tof_bin_index, num_tof_bins, lor_dependent_sigma_tof, lor_dependent_tofcenter_offset)

%--------------------------------------------------------------------------
% forward projection 
%--------------------------------------------------------------------------

image_size_xy = img_dim(1);
image_size_z = img_dim(3);

image_recon = reshape(image_recon, [image_size_xy, image_size_xy, image_size_z]);
image_recon = permute(image_recon, [3, 2, 1]);

image_recon = single(image_recon(:));

sino_init=zeros(size(sino_init));
sino_init = single(sino_init(:));

device_id = int32(0);
threads = int32(64);

% ✅ 类型安全转换，与 C 头文件签名严格对应
tof_bin_index              = int16(tof_bin_index);           % const short*
num_tof_bins               = int16(num_tof_bins);            % short
lor_dependent_sigma_tof    = uint8(lor_dependent_sigma_tof); % unsigned char
lor_dependent_tofcenter_offset = uint8(lor_dependent_tofcenter_offset); % unsigned char
nlors                      = uint64(nlors);                  % size_t
device_id                  = int32(device_id);
threads                    = int32(threads);

% 严格请求 10 个输出参数，对应 C 头文件中的 10 个指针
[~, ~, ~, ~, ~, sino_fwd, ~, ~, ~, ~] = calllib(libname, 'joseph3d_tof_lm_fwd', xstart, xend, image_recon, img_origin, voxel_size, sino_init, nlors, img_dim, tofbin_width, sigma_tof, tofcenter_offset, n_sigmas_TOF, tof_bin_index, num_tof_bins, lor_dependent_sigma_tof, lor_dependent_tofcenter_offset, device_id, threads);

sino_fwd = double(sino_fwd(:));

end