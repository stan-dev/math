#ifndef STAN_MATH_OPENCL_KERNELS_ADD_HPP
#define STAN_MATH_OPENCL_KERNELS_ADD_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static const char *add_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Matrix addition on the OpenCL device
     *
     * @param[out] C Output matrix.
     * @param[in] A LHS of matrix addition.
     * @param[in] B RHS of matrix addition.
     * @param rows Number of rows for matrix A.
     * @param cols Number of cols for matrix A.
     * @note Code is a <code>const char*</code> held in
     * <code>add_kernel_code.</code>
     * This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void add(__global double *C, __global double *A,
                      __global double *B, unsigned int rows,
                      unsigned int cols) {
      const int i = get_global_id(0);
      const int j = get_global_id(1);
      if (i < rows && j < cols) {
        C(i, j) = A(i, j) + B(i, j);
      }
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/add.hpp add() \endlink
 */
const global_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, int, int> add(
    "add", add_kernel_code);
// \cond
static const char *add_batch_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Adds all the matrices in a batch
     *
     * @param[in,out] A buffer containing the entire batch.
     * @param rows Number of rows for a single matrix in the batch.
     * @param cols Number of cols for a single matrix in the batch.
     * @param batch_size Number of matrices in the batch.
     * @note Code is a <code>const char*</code> held in
     * <code>add_batch_kernel_code.</code>
     * This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void add_batch(__global double *A, unsigned int rows,
                      unsigned int cols, unsigned int batch_size) {
      const int i = get_global_id(0);
      const int j = get_global_id(1);
      if (i < rows && j < cols) {
        for(int k = 1;k < batch_size;k++){
          A_batch(i, j, 0) += A_batch(i, j, k);
        }
      }
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/add.hpp add_batch() \endlink
 */
const global_range_kernel<cl::Buffer, int, int, int> add_batch(
    "add_batch", add_batch_kernel_code);
}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
