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
     * @param[in, out] V A temporary matrix for storing partial calculations.
     * @param rows The number of rows for A.
     * @note Code is a <code>const char*</code> held in
     * <code>lower_tri_inverse_step1_kernel_code.</code>
     *  Used in math/gpu/lower_tri_inverse.hpp.
     *  This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void lower_tri_inverse_step1(__global double* A,
                                          __global double* V, int rows) {
      int index = get_local_id(0);
      int group = get_group_id(0);
      int block_size = get_local_size(0);
      int offset = group * block_size;
      for (int j = 0; j < block_size; j++) {
        if (index == j) {
          V[group * block_size * block_size + j * block_size + index] = 1.0;
        } else {
          V[group * block_size * block_size + j * block_size + index] = 0.0;
        }
      }
      barrier(CLK_LOCAL_MEM_FENCE);
      for (int k = 0; k < block_size; k++) {
        double factor = A(offset + k, offset + k);
        if (index <= k) {
          V[group * block_size * block_size + index * block_size + k] /= factor;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int i = k + 1; i < block_size; i++) {
          factor = A(offset + i, offset + k);
          if (index < i) {
            V[group * block_size * block_size + index * block_size + i]
                -= V[group * block_size * block_size + index * block_size + k]
                   * factor;
          }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      for (int j = 0; j < block_size; j++) {
        A(offset + index, offset + j)
            = V[group * block_size * block_size + j * block_size + index];
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
