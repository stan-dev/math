#ifndef STAN_MATH_GPU_KERNELS_IDENTITY_HPP
#define STAN_MATH_GPU_KERNELS_IDENTITY_HPP
#ifdef STAN_OPENCL

#include <stan/math/gpu/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
const char* identity_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Makes an identity matrix on the GPU
     *
     * @param[in,out] A The identity matrix output.
     * @param rows The number of rows for A.
     * @param cols The number of cols for A.
     * @note Code is a <code>const char*</code> held in
     * <code>identity_kernel_code.</code>
     *  Used in math/gpu/identity_opencl.hpp.
     *  This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void identity(__global double* A, unsigned int rows,
                           unsigned int cols) {
      int i = get_global_id(0);
      int j = get_global_id(1);
      if (i < rows && j < cols) {
        if (i == j) {
          A(i, j) = 1.0;
        } else {
          A(i, j) = 0.0;
        }
      }
    }
    // \cond
);
// \endcond
// \cond
const char* batch_identity_kernel_code = STRINGIFY(
    // \endcond

    /**
     * Makes a batch of smaller identity matrices inside the input matrix
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
      int batch_id = get_global_id(0);
      int batch_row = get_global_id(1);
      int batch_col = get_global_id(2);
      int index = batch_id * batch_rows * batch_rows + batch_col * batch_rows
                  + batch_row;
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

/**
 * See the docs for \link kernels/identity.hpp identity() \endlink
 */
const global_range_kernel<cl::Buffer, int, int> identity("identity",
                                                         identity_kernel_code);

/**
 * See the docs for \link kernels/identity.hpp batch_identity() \endlink
 */
const global_range_kernel<cl::Buffer, int, int> batch_identity(
    "batch_identity", batch_identity_kernel_code);

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
