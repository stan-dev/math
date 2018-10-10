#ifndef STAN_MATH_GPU_KERNELS_LOWER_TRI_INVERSE_STEP1_HPP
#define STAN_MATH_GPU_KERNELS_LOWER_TRI_INVERSE_STEP1_HPP
#ifdef STAN_OPENCL

#include <stan/math/gpu/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
const char* lower_tri_inverse_step1_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Calculates the inverse with no blocking
     *
     * @param[in,out] A The input matrix.
     * @param[in, out] tmp_inv A temporary matrix for storing partial calculations.
     * @param rows The number of rows for A.
     * @note Code is a <code>const char*</code> held in
     * <code>lower_tri_inverse_step1_kernel_code.</code>
     *  Used in math/gpu/lower_tri_inverse.hpp.
     *  This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void lower_tri_inverse_step1(__global double* A,
                                          __global double* tmp_inv, int rows) {
      int index = get_local_id(0);
      int group = get_group_id(0);
      int block_size = get_local_size(0);
      int offset = group * block_size;
      int tmp_offset = group * block_size * block_size + index * block_size;
      for (int k = 0; k < block_size; k++) {
        double factor = A(offset + k, offset + k);
        if (index <= k) {
          tmp_inv[tmp_offset + k] /= factor;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int i = max(k + 1, index); i < block_size; i++) {
          factor = A(offset + i, offset + k);
            tmp_inv[tmp_offset + i]
                -= tmp_inv[tmp_offset + k]
                   * factor;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      for (int j = 0; j < block_size; j++) {
        A(offset + j, offset + index)
            = tmp_inv[tmp_offset + j];
      }
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/lower_tri_inverse_step1.hpp add() \endlink
 */
const local_range_kernel<cl::Buffer, cl::Buffer, int> lower_tri_inverse_step1(
    "lower_tri_inverse_step1", lower_tri_inverse_step1_kernel_code,
    {{"THREAD_BLOCK_SIZE", 32}});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
