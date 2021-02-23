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
static const std::string rep_matrix_rev_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
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
    __kernel void rep_matrix_rev(__global double* A_adj, __global double* B_adj,
                             unsigned int B_rows, unsigned int B_cols,
                             unsigned int view_B) {
      const int gid_i = get_global_id(0);
      const int gid_j = get_global_id(1);
      const int gsize_i = get_global_size(0);
      const int gsize_j = get_global_size(1);
      double tmp = 0;
      //      j_start = contains_nonzero(view_B, LOWER) ? gid_j :
      for (int j = gid_j; j < B_cols; j += gsize_j) {
        int i_start
            = contains_nonzero(view_B, UPPER)
                  ? gid_i
                  : ((j - gid_i + gsize_i - 1) / gsize_i) * gsize_i + gid_i;
        int i_end = contains_nonzero(view_B, LOWER) ? B_rows : j+1;
        for (int i = i_start; i < i_end; i += gsize_i) {
          tmp += B_adj[j * B_rows + i];
        }
      }
      A_adj[gid_j * gsize_i + gid_i] += tmp;

//      if (i < A_rows && j < A_cols) {
//        double val = 0;
//        if (B_cols == 1 && B_rows == 1) {
//          val = B[0];
//        } else if (B_cols == 1) {
//          val = B[i];
//        } else if (B_rows == 1) {
//          val = B[j];
//        }
//        if ((contains_nonzero(view_A, LOWER) && j <= i)
//            || (contains_nonzero(view_A, UPPER) && j >= i)) {
//          A[j * A_rows + i] = val;
//        }
//      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/rep_matrix.hpp rep_matrix_rev() \endlink
 */
const kernel_cl<in_out_buffer, in_buffer, int, int, matrix_cl_view>
    rep_matrix_rev("rep_matrix_rev",
               {view_kernel_helpers, rep_matrix_rev_kernel_code});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
