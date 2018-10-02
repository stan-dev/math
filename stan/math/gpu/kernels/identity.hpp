#ifndef STAN_MATH_GPU_KERNELS_IDENTITY_HPP
#define STAN_MATH_GPU_KERNELS_IDENTITY_HPP
#ifdef STAN_OPENCL

#include <stan/math/gpu/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
const char *identity_kernel_code = STRINGIFY(
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
    __kernel void identity(__global double *A, unsigned int rows,
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

/**
 * See the docs for \link kernels/transpose.hpp transpose() \endlink
 */
const global_range_kernel<cl::Buffer, int, int> identity("identity",
                                                         identity_kernel_code);

} // namespace opencl_kernels
} // namespace math
} // namespace stan
#endif
#endif
