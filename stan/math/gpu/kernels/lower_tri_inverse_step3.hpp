#ifndef STAN_MATH_GPU_KERNELS_LOWER_TRI_INVERSE_STEP3_HPP
#define STAN_MATH_GPU_KERNELS_LOWER_TRI_INVERSE_STEP3_HPP
#ifdef STAN_OPENCL

#include <stan/math/gpu/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
const char* lower_tri_inverse_step3_kernel_code = STRINGIFY(
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
__kernel void lower_tri_inverse_step3(
        __global double* A,
        __global int* sizes,
        __global double* temp,
        int repeat,
        int remainder,
        int part_size_fixed,
        int rows) {
 
}
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/matrix_multiply.hpp add() \endlink
 */
const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, int, int, int, int>
    lower_tri_inverse_step3("lower_tri_inverse_step3", lower_tri_inverse_step3_kernel_code);

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
