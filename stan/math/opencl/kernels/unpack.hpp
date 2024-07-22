#ifndef STAN_MATH_OPENCL_KERNELS_UNPACK_HPP
#define STAN_MATH_OPENCL_KERNELS_UNPACK_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static constexpr const char* unpack_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     * Unpacks a packed triangular matrix to a flat
     * matrix
     *
     * @param[out] B flat matrix
     * @param[in] A packed buffer
     * @param rows number of columns for matrix B
     * @param cols number of columns for matrix B
     * @param view parameter that defines the triangularity of the
     * input matrix
     *  LOWER - lower triangular
     *  UPPER - upper triangular
     * if the view parameter is not specified
     * @note Code is a <code>const char*</code> held in
     * <code>unpack_kernel_code.</code>
     * This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void unpack(__global double* B, __global double* A,
                         unsigned int rows, unsigned int cols,
                         unsigned int view) {
      int i = get_global_id(0);
      int j = get_global_id(1);
      if (i < rows && j < cols) {
        // the packed matrices are stored in row major
        if (view == LOWER) {
          const int column_offset = j * rows - (j * (j - 1)) / 2;
          const int row_offset = (i - j);
          if (j <= i) {
            B(i, j) = A[column_offset + row_offset];
          } else {
            B(i, j) = 0.0;
          }
        } else {
          const int column_offset = j * (j + 1) / 2;
          if (j >= i) {
            B(i, j) = A[column_offset + i];
          } else {
            B(i, j) = 0.0;
          }
        }
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/unpack.hpp unpack() \endlink
 */
const kernel_cl<out_buffer, in_buffer, int, int, matrix_cl_view> unpack(
    "unpack", {indexing_helpers, unpack_kernel_code});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
