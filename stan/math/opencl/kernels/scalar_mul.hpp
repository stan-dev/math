#ifndef STAN_MATH_OPENCL_KERNELS_SCALAR_MUL_HPP
#define STAN_MATH_OPENCL_KERNELS_SCALAR_MUL_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static const std::string scalar_mul_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Multiplication of the matrix A with a scalar
     *
     * @param[out] A output matrix
     * @param[in] B input matrix
     * @param[in] scalar the value with which to multiply A
     * @param[in] rows the number of rows in A
     * @param[in] cols the number of columns in A
     * @param[in] view triangular view of the input matrix to use
     */
    __kernel void scalar_mul(__global double *A, const __global double *B,
                             const double scalar, const unsigned int rows,
                             const unsigned int cols, int view) {
      int i = get_global_id(0);
      int j = get_global_id(1);
      if (i < rows && j < cols) {
        if (!((!contains_nonzero(view, LOWER) && j < i)
              || (!contains_nonzero(view, UPPER) && j > i))) {
          A(i, j) = B(i, j) * scalar;
        } else {
          A(i, j) = 0;
        }
      }
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/scalar_mul.hpp add() \endlink
 */
const kernel_cl<out_buffer, in_buffer, double, int, int, matrix_cl_view>
    scalar_mul("scalar_mul",
               {indexing_helpers, view_kernel_helpers, scalar_mul_kernel_code});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
