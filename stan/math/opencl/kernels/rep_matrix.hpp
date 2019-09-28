#ifndef STAN_MATH_OPENCL_KERNELS_REP_MATRIX_HPP
#define STAN_MATH_OPENCL_KERNELS_REP_MATRIX_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static const std::string rep_matrix_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Creates a matrix from a matrix_cl of size 1x1 by
     * replicating the single value or by replicating the
     * vector or row_vector input.
     *
     * @param[out] A result matrix
     * @param[in] B input matrix (1x1, vector or row_vector)
     * @param A_rows Number of rows for matrix A
     * @param A_cols Number of columns for matrix A
     * @param B_rows Number of rows for matrix B
     * @param B_cols Number of columns for matrix B
     * @param view_A triangular part of matrix A to use
     *
     * @note Code is a string held in <code>rep_matrix_kernel_code.</code>
     * This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void rep_matrix(__global double* A, __global double* B,
                             unsigned int A_rows, unsigned int A_cols,
                             unsigned int B_rows, unsigned int B_cols,
                             unsigned int view_A) {
      const int i = get_global_id(0);
      const int j = get_global_id(1);
      if (i < A_rows && j < A_cols) {
        double val = 0;
        if (B_cols == 1 && B_rows == 1) {
          val = B[0];
        } else if (B_cols == 1) {
          val = B[i];
        } else if (B_rows == 1) {
          val = B[j];
        }
        if ((contains_nonzero(view_A, LOWER) && j <= i)
            || (contains_nonzero(view_A, UPPER) && j >= i)) {
          A[j * A_rows + i] = val;
        }
      }
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/rep_matrix.hpp rep_matrix() \endlink
 */
const kernel_cl<out_buffer, in_buffer, int, int, int, int, matrix_cl_view>
    rep_matrix("rep_matrix",
               {indexing_helpers, view_kernel_helpers, rep_matrix_kernel_code});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
