#ifndef STAN_MATH_OPENCL_KERNELS_FILL_STRICT_TRI_HPP
#define STAN_MATH_OPENCL_KERNELS_FILL_STRICT_TRI_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static const std::string fill_strict_tri_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     * Stores constant in the triangular part of a matrix
     * on the OpenCL device. Supports writing constants
     * to the lower and upper triangular. The input
     * matrix is unchanged if used with any other view.
     *
     * @param[out] A matrix
     * @param val value to replicate in the matrix
     * @param rows Number of rows for matrix A
     * @param cols Number of columns for matrix A
     * @param view_A triangular part of matrix A to use
     *
     * @note Code is a <code>const char*</code> held in
     * <code>fill_strict_tri_kernel_code.</code>
     * This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void fill_strict_tri(__global double* A, double val,
                                  unsigned int rows, unsigned int cols,
                                  unsigned int view_A) {
      const int i = get_global_id(0);
      const int j = get_global_id(1);
      if (i < rows && j < cols) {
        if ((view_A == LOWER && j < i) || (view_A == UPPER && j > i)
            || view_A == ENTIRE) {
          A(i, j) = val;
        }
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/fill_strict_tri.hpp
 * fill_strict_tri_kernel_code() \endlink
 */
const kernel_cl<out_buffer, double, int, int, matrix_cl_view> fill_strict_tri(
    "fill_strict_tri",
    {indexing_helpers, view_kernel_helpers, fill_strict_tri_kernel_code});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
