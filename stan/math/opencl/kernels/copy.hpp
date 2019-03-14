#ifndef STAN_MATH_OPENCL_KERNELS_COPY_HPP
#define STAN_MATH_OPENCL_KERNELS_COPY_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <algorithm>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static const char *copy_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Copy one matrix to another
     * @param[in] A The matrix to copy.
     * @param[out] B The matrix to copy A to.
     * @param rows The number of rows in A.
     * @param cols The number of cols in A.
     * @note Code is a <code>const char*</code> held in
     * <code>copy_kernel_code.</code>
     * Kernel used in math/opencl/matrix_cl.hpp.
     *  This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void copy(__global double *A, __global double *B,
                       unsigned int rows, unsigned int cols,
                       triangular_view tri_view) {
      int i = get_global_id(0);
      int j = get_global_id(1);
      if (i < rows && j < cols) {
        if (tri_view == LOWER && j <= i) {
          B(i, j) = A(i, j);
        } else if (tri_view == LOWER) {
          B(i, j) = 0;
        } else if (tri_view == UPPER && j >= i) {
          B(i, j) = A(i, j);
        } else if (tri_view == UPPER) {
          B(i, j) = 0;
        } else {
          // tri_view == ENTIRE
          B(i, j) = A(i, j);
        }
      }
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/copy.hpp copy() \endlink
 */
const global_range_kernel<cl::Buffer, cl::Buffer, int, int,
                                     triangular_view_CL> copy(
    "copy", copy_kernel_code);

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
