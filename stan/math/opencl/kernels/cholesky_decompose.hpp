#ifndef STAN_MATH_OPENCL_KERNELS_CHOLESKY_DECOMPOSE_HPP
#define STAN_MATH_OPENCL_KERNELS_CHOLESKY_DECOMPOSE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static constexpr const char* cholesky_decompose_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     * Calculates the Cholesky Decomposition of a matrix on an OpenCL
     *
     * This kernel is run with threads organized in one dimension and
     * in a single thread block. The kernel is best suited for
     * small input matrices as it only utilizes a single streaming
     * multiprocessor. The kernels is used as a part of a blocked
     * cholesky decompose.
     *
     * @param[in, out] A The input matrix and the result of the cholesky
     *  decomposition
     * @param rows The number of rows for A and B.
     * @param row_j Local memory buffer of rows doubles that holds row j.
     * @note Code is a <code>const char*</code> held in
     * <code>cholesky_decompose_kernel_code.</code>
     *  Used in math/opencl/cholesky_decompose.hpp.
     *  This kernel uses the helper macros available in helpers.cl.
     *
     */
    __kernel void cholesky_decompose(__global double* A, int rows,
                                     __local double* row_j) {
      const int local_index = get_local_id(0);
      // The following code is the sequential version of the inplace
      // cholesky decomposition. Only the innermost loops are parallelized. The
      // rows are processed sequentially. This loop process all the rows:
      for (int j = 0; j < rows; j++) {
        if (local_index == j) {
          double sum = 0.0;
          for (int k = 0; k < j; k++) {
            const double a = A(j, k);
            row_j[k] = a;
            sum = sum + a * a;
          }
          A(j, j) = sqrt(A(j, j) - sum);
          row_j[j] = A(j, j);
        }
        // Each work item accesses only its own row of A and reads row j from
        // local memory, so the barriers need only a local memory fence.
        barrier(CLK_LOCAL_MEM_FENCE);
        if (local_index < j) {
          A(local_index, j) = 0.0;
        } else if (local_index > j) {
          double sum = 0.0;
          for (int k = 0; k < j; k++)
            sum = sum + A(local_index, k) * row_j[k];
          A(local_index, j) = (A(local_index, j) - sum) / row_j[j];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/cholesky_decompose.hpp cholesky_decompose()
 * \endlink
 */
const kernel_cl<in_out_buffer, int, cl::LocalSpaceArg> cholesky_decompose(
    "cholesky_decompose", {indexing_helpers, cholesky_decompose_kernel_code});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
