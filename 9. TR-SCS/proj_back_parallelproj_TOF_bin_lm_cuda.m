function Ratio_factor_image = proj_back_parallelproj_TOF_bin_lm_cuda(libname, xstart, xend, image_recon, img_origin, voxel_size, Ratio_factor_sino, nlors, img_dim, tofbin_width, sigma_tof, tofcenter_offset, n_sigmas_TOF, tof_bin_index, num_tof_bins, lor_dependent_sigma_tof, lor_dependent_tofcenter_offset)

%--------------------------------------------------------------------------
% back projection
%--------------------------------------------------------------------------

image_size_xy = img_dim(1);
image_size_z = img_dim(3);

image_recon = zeros(size(image_recon));

image_recon = single(image_recon(:));
Ratio_factor_sino = single(Ratio_factor_sino(:));

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

[~, ~, Ratio_factor_image, ~, ~, ~, ~, ~, ~, ~] = calllib(libname, 'joseph3d_tof_lm_back', xstart, xend, image_recon, img_origin, voxel_size, Ratio_factor_sino, nlors, img_dim, tofbin_width, sigma_tof, tofcenter_offset, n_sigmas_TOF, tof_bin_index, num_tof_bins, lor_dependent_sigma_tof, lor_dependent_tofcenter_offset, device_id, threads);

Ratio_factor_image = reshape(Ratio_factor_image, [image_size_z, image_size_xy, image_size_xy]);
Ratio_factor_image = permute(Ratio_factor_image, [3, 2, 1]);

Ratio_factor_image = double(Ratio_factor_image(:));
end