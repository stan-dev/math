#ifndef STAN_MATH_OPENCL_KERNELS_CONSTANTS_HPP
#define STAN_MATH_OPENCL_KERNELS_CONSTANTS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static const std::string constants_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Stores constant in the matrix on the OpenCL device.
     * Supports writing constants to the lower and upper triangular or
     * the whole matrix.
     *
     * @param[out] A matrix
     * @param val value to replicate in the matrix
     * @param rows Number of rows for matrix A
     * @param cols Number of columns for matrix A
     * @param view_A triangular part of matrix A to use
     *
     * @note Code is a <code>const char*</code> held in
     * <code>constant_kernel_code.</code>
     * This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void constants(__global double* A, double val, unsigned int rows,
                            unsigned int cols, unsigned int view_A) {
      const int i = get_global_id(0);
      const int j = get_global_id(1);
      if (i < rows && j < cols) {
        if ((contains_nonzero(view_A, LOWER) && j <= i)
            || (contains_nonzero(view_A, UPPER) && j >= i)) {
          A(i, j) = val;
        }
      }
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/constants.hpp constants() \endlink
 */
const kernel_cl<out_buffer, double, int, int, matrix_cl_view> constants(
    "constants",
    {indexing_helpers, view_kernel_helpers, constants_kernel_code});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
