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
        __global double* V,
        int rows) {
      int index = get_local_id(0);
      int group = get_group_id(0);
      int block_size = get_local_size(0);
      int offset = group*block_size;
      //__local double V[THREAD_BLOCK_SIZE][THREAD_BLOCK_SIZE];      
      int i = index;
      for(int j=0;j<block_size;j++){
        if(i==j){
            V[group*block_size*block_size+j*block_size + i] = 1.0;
        }else{
            V[group*block_size*block_size+j*block_size + i] = 0.0;
        }
      }
      barrier(CLK_LOCAL_MEM_FENCE);
      for(int k=0;k<block_size;k++){
        
        double factor = A(offset+k,offset+k); 
        int j = index;
        if (j <= k) {
            V[group*block_size*block_size+j*block_size + k] /= factor;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        for(int i=k+1;i<block_size;i++){
            factor = A(offset+i,offset+k);
            j = index;
            if (j < i) {
                V[group*block_size*block_size+j*block_size + i] -= V[group*block_size*block_size+j*block_size + k]*factor;
            }
        }        
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      barrier(CLK_LOCAL_MEM_FENCE);
      i = index;
      for(int j=0;j<block_size;j++){
        A(offset+i,offset+j) = V[group*block_size*block_size+j*block_size + i];
      }
      
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/matrix_multiply.hpp add() \endlink
 */
const local_range_kernel<cl::Buffer, cl::Buffer, int>
    lower_tri_inverse_step1("lower_tri_inverse_step1", lower_tri_inverse_step1_kernel_code,
    {{"THREAD_BLOCK_SIZE", 32}});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
