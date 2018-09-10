#ifndef STAN_MATH_GPU_KERNELS_LOWER_TRI_INVERSE_STEP1_HPP
#define STAN_MATH_GPU_KERNELS_LOWER_TRI_INVERSE_STEP1_HPP
#ifdef STAN_OPENCL

#include <stan/math/gpu/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
const char* lower_tri_inverse_step1_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Calculates 
     *
     * @param[in,out] A The identity matrix output.
     * @param rows The number of rows for A.
     * @param cols The number of cols for A.
     * @note Code is a <code>const char*</code> held in
     * <code>identity_kernel_code.</code>
     *  Used in math/gpu/identity_opencl.hpp.
     *  This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void lower_tri_inverse_step1(
        __global double* A,
        int rows) {
      int i = get_global_id(0);
      double a = 0;
      double b = 0;
      double c = 0;
      if (2*i < rows) {
          a = 1.0/A(2*i, 2*i);
          A(2*i, 2*i) = a;
      }
      if (2*i+1 < rows) {
          b = 1.0/A(2*i+1, 2*i+1);
          c = A(2*i+1, 2*i);
          A(2*i+1, 2*i+1) = b;
          A(2*i+1, 2*i) = -a*b*c;
      }
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/matrix_multiply.hpp add() \endlink
 */
const global_range_kernel<cl::Buffer, int>
    lower_tri_inverse_step1("lower_tri_inverse_step1", lower_tri_inverse_step1_kernel_code);

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
