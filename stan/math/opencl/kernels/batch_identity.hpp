#ifndef STAN_MATH_OPENCL_KERNELS_IDENTITY_HPP
#define STAN_MATH_OPENCL_KERNELS_IDENTITY_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static constexpr const char* batch_identity_kernel_code = STRINGIFY(
    // \endcond

    /** \ingroup opencl_kernels
     * Makes a batch of smaller identity matrices inside the input matrix
     *
     * This kernel operates inplace on the matrix A, filling it with smaller
     *  identity matrices with a size of batch_rows x batch_rows.
     *  This kernel expects a 3D organization of threads:
     *      1st dim: the number of matrices in the batch.
     *      2nd dim: the number of cols/rows in batch matrices.
     *      3rd dim: the number of cols/rows in batch matrices.
     *  Each thread in the organization assigns a single value in the batch.
     *  In order to create a batch of 3 matrices the size of NxN you need
     *  to run the kernel batch_identity(A, N, 3*N*N) with (3, N, N) threads.
     *  The special case of batch_identity(A, N, N*N) executed on
     *  (1, N, N) threads creates a single identity matrix the size of NxN and
     *  is therefore equal to the basic identity kernel.
     *
     * @param[in,out] A The batched identity matrix output.
     * @param batch_rows The number of rows/cols for the smaller matrices in the
     * batch
     * @param size The size of A.
     * @note Code is a <code>const char*</code> held in
     * <code>identity_kernel_code.</code>
     *  This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void batch_identity(__global double* A, unsigned int batch_rows,
                                 unsigned int size) {
      // The ID of the matrix in the batch the thread is assigned to
      int batch_id = get_global_id(0);
      // The row and column of the matrix in the batch
      int batch_row = get_global_id(1);
      int batch_col = get_global_id(2);
      int index = batch_id * batch_rows * batch_rows + batch_col * batch_rows
                  + batch_row;
      // Check for potential overflows of A.
      if (index < size) {
        if (batch_row == batch_col) {
          A[index] = 1.0;
        } else {
          A[index] = 0.0;
        }
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/batch_identity.hpp batch_identity() \endlink
 */
const kernel_cl<out_buffer, int, int> batch_identity(
    "batch_identity", {indexing_helpers, batch_identity_kernel_code});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
