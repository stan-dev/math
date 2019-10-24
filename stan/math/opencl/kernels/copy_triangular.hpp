#ifndef STAN_MATH_OPENCL_KERNELS_COPY_TRIANGULAR_HPP
#define STAN_MATH_OPENCL_KERNELS_COPY_TRIANGULAR_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static const std::string copy_triangular_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Copies the lower or upper
     * triangular of the source matrix to
     * the destination matrix.
     * Both matrices are stored on the OpenCL device.
     *
     * @param[out] A Output matrix to copy triangular to.
     * @param[in] B The matrix to copy the triangular from.
     * @param rows The number of rows of B.
     * @param cols The number of cols of B.
     * @param view determines
     * which part of the matrix to copy:
     *  ENTIRE: copies entire matrix
     *  LOWER: copies the lower triangular
     *  UPPER: copies the upper triangular
     *  DIAGONAL: copies the diagonal
     * @note Code is a <code>const char*</code> held in
     * <code>copy_triangular_kernel_code.</code>
     * Used in math/opencl/copy_triangular_opencl.hpp.
     *  This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void copy_triangular(__global double *A, __global double *B,
                                  unsigned int rows, unsigned int cols,
                                  unsigned int view) {
      int i = get_global_id(0);
      int j = get_global_id(1);
      if (i < rows && j < cols) {
        if ((contains_nonzero(view, LOWER) && j <= i)
            || (contains_nonzero(view, UPPER) && j >= i) || j == i) {
          A(i, j) = B(i, j);
        } else {
          A(i, j) = 0;
        }
      }
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/copy_triangular.hpp copy_triangular() \endlink
 */
const kernel_cl<out_buffer, in_buffer, int, int, matrix_cl_view>
    copy_triangular("copy_triangular", {indexing_helpers, view_kernel_helpers,
                                        copy_triangular_kernel_code});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
