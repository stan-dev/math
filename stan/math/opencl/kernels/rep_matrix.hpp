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
     * Implements reverse pass of rep_matrix.
     *
     * @param[in,out] A_adj adjoint of argument matrix
     * @param[in] B_adj adjoint of result matrix
     * @param B_rows Number of rows for matrix B
     * @param B_cols Number of columns for matrix B
     * @param view_B triangular part of matrix B to use
     */
    __kernel void rep_matrix_rev(__global double* A_adj, __global double* B_adj,
                                 unsigned int B_rows, unsigned int B_cols,
                                 unsigned int view_B) {
      const int gid_i = get_global_id(0);
      const int gid_j = get_global_id(1);
      const int gsize_i = get_global_size(0);
      const int gsize_j = get_global_size(1);
      double tmp = 0;
      for (int j = gid_j; j < B_cols; j += gsize_j) {
        int i_start
            = contains_nonzero(view_B, UPPER)
                  ? gid_i
                  : ((j - gid_i + gsize_i - 1) / gsize_i) * gsize_i + gid_i;
        int i_end = contains_nonzero(view_B, LOWER) ? B_rows : j + 1;
        for (int i = i_start; i < i_end; i += gsize_i) {
          tmp += B_adj[j * B_rows + i];
        }
      }
      A_adj[gid_j * gsize_i + gid_i] += tmp;
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
