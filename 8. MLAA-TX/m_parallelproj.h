/* m_parallelproj.h - MATLAB compatible C wrapper for parallelproj */

#ifndef M_PARALLELPROJ_H
#define M_PARALLELPROJ_H

#include <stddef.h> /* 替换 <cstddef> */

/* 屏蔽原始的导出宏，防止 MATLAB 找不到 parallelproj_export.h */
#define PARALLELPROJ_API 

#ifdef __cplusplus
extern "C" {
#endif

  PARALLELPROJ_API int parallelproj_cuda_enabled(void);
  PARALLELPROJ_API const char *parallelproj_version(void);
  PARALLELPROJ_API int parallelproj_version_major(void);
  PARALLELPROJ_API int parallelproj_version_minor(void);
  PARALLELPROJ_API int parallelproj_version_patch(void);

  /* 去除了 std:: 命名空间，去除了 = 0 和 = 64 这种 C++ 默认参数 */
  PARALLELPROJ_API void joseph3d_fwd(const float *lor_start, const float *lor_end, const float *image, const float *image_origin, const float *voxel_size, float *projection_values, size_t num_lors, const int *image_dim, int device_id, int threads_per_block);

  PARALLELPROJ_API void joseph3d_back(const float *lor_start, const float *lor_end, float *image, const float *image_origin, const float *voxel_size, const float *projection_values, size_t num_lors, const int *image_dim, int device_id, int threads_per_block);

  PARALLELPROJ_API void joseph3d_tof_sino_fwd(const float *lor_start, const float *lor_end, const float *image, const float *image_origin, const float *voxel_size, float *projection_values, size_t num_lors, const int *image_dim, float tof_bin_width, const float *tof_sigma, const float *tof_center_offset, float num_sigmas, short num_tof_bins, unsigned char is_lor_dependent_tof_sigma, unsigned char is_lor_dependent_tof_center_offset, int device_id, int threads_per_block);

  PARALLELPROJ_API void joseph3d_tof_sino_back(const float *lor_start, const float *lor_end, float *image, const float *image_origin, const float *voxel_size, const float *projection_values, size_t num_lors, const int *image_dim, float tof_bin_width, const float *tof_sigma, const float *tof_center_offset, float num_sigmas, short num_tof_bins, unsigned char is_lor_dependent_tof_sigma, unsigned char is_lor_dependent_tof_center_offset, int device_id, int threads_per_block);

  PARALLELPROJ_API void joseph3d_tof_lm_fwd(const float *lor_start, const float *lor_end, const float *image, const float *image_origin, const float *voxel_size, float *projection_values, size_t num_events, const int *image_dim, float tof_bin_width, const float *tof_sigma, const float *tof_center_offset, float num_sigmas, const short *tof_bin_index, short num_tof_bins, unsigned char is_lor_dependent_tof_sigma, unsigned char is_lor_dependent_tof_center_offset, int device_id, int threads_per_block);

  PARALLELPROJ_API void joseph3d_tof_lm_back(const float *lor_start, const float *lor_end, float *image, const float *image_origin, const float *voxel_size, const float *projection_values, size_t num_events, const int *image_dim, float tof_bin_width, const float *tof_sigma, const float *tof_center_offset, float num_sigmas, const short *tof_bin_index, short num_tof_bins, unsigned char is_lor_dependent_tof_sigma, unsigned char is_lor_dependent_tof_center_offset, int device_id, int threads_per_block);

#ifdef __cplusplus
}
#endif

#endif
