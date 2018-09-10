#ifndef STAN_MATH_GPU_KERNELS_LOWER_TRI_INVERSE_HPP
#define STAN_MATH_GPU_KERNELS_LOWER_TRI_INVERSE_HPP
#ifdef STAN_OPENCL

#include <stan/math/gpu/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
const char* lower_tri_inverse_step2_small_kernel_code = STRINGIFY(
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
    __kernel void lower_tri_inverse_step2_small(
        __global double* A,
        int temp_rows,
        int rows) {
            int i = get_global_id(0); 
            int j = get_global_id(1);
            int t = get_global_id(2);
            double acc = 0.0;
            int offset = t*temp_rows*2;
            __local double temp[THREAD_BLOCK_SIZE][THREAD_BLOCK_SIZE];
            for(int k = 0; k < temp_rows; k++){
                acc += A(offset+temp_rows+i,offset+temp_rows+k)*A(offset+temp_rows+k,offset+j);
            }
            temp[i][j]=acc;
            barrier(CLK_LOCAL_MEM_FENCE);
            acc = 0.0;
            for(int k = 0; k < temp_rows; k++){
                acc += temp[i][k]*A(offset+k, offset+j);
            }
            A(offset+temp_rows+i, offset+j) = -acc;

    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/matrix_multiply.hpp add() \endlink
 */
const local_range_kernel<cl::Buffer, int, int>
    lower_tri_inverse_step2_small("lower_tri_inverse_step2_small", lower_tri_inverse_step2_small_kernel_code,
    {{"THREAD_BLOCK_SIZE", 32}});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
