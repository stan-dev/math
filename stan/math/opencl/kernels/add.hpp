#ifndef STAN_MATH_OPENCL_KERNELS_ADD_HPP
#define STAN_MATH_OPENCL_KERNELS_ADD_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static constexpr const char *add_batch_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     * Sums a batch of matrices. Buffer A contains
     * batch_size matrices of size rows x cols. All elements
     * at matching indices are summed up and stored to the
     * resulting matrix B.
     *
     * @param[out] B buffer of the result matrix.
     * @param[in] A buffer containing the entire batch.
     * @param rows Number of rows for a single matrix in the batch.
     * @param cols Number of cols for a single matrix in the batch.
     * @param batch_size Number of matrices in the batch.
     * @note Code is a <code>const char*</code> held in
     * <code>add_batch_kernel_code.</code>
     * This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void add_batch(__global double *B, __global double *A,
                            unsigned int rows, unsigned int cols,
                            unsigned int batch_size) {
      const int i = get_global_id(0);
      const int j = get_global_id(1);
      if (i < rows && j < cols) {
        double temp = 0.0;
        for (int k = 0; k < batch_size; k++) {
          temp += A_batch(i, j, k);
        }
        B(i, j) = temp;
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/add.hpp add_batch() \endlink
 */
const kernel_cl<out_buffer, in_buffer, int, int, int> add_batch(
    "add_batch", {indexing_helpers, add_batch_kernel_code});
}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
