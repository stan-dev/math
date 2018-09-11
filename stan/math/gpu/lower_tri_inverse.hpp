#ifndef STAN_MATH_GPU_LOWER_TRI_INVERSE_HPP
#define STAN_MATH_GPU_LOWER_TRI_INVERSE_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/matrix_gpu.hpp>
#include <stan/math/gpu/kernels/lower_tri_inverse_step1.hpp>
#include <stan/math/gpu/kernels/lower_tri_inverse_step2.hpp>
#include <stan/math/gpu/kernels/lower_tri_inverse_step3.hpp>
#include <stan/math/gpu/err/check_square.hpp>
#include <string>
#include <vector>

namespace stan {
  namespace math {
    /**
     * Computes the inverse of the lower triangular matrix 
     * that resides in the GPU global memory
     * 
     * @param A matrix on the GPU 
     *
     * @return the inverse of A
     * 
     * @throw <code>std::invalid_argument</code> if the matrix
     *    is not square
     */
    inline matrix_gpu lower_triangular_inverse(const matrix_gpu & A) {
      check_square("lower_triangular_inverse (GPU)", "A", A);
      matrix_gpu inv(A);
      try {
        matrix_gpu temp((inv.rows()+1)/2,  (inv.cols()+1)/2);
        int thread_block_size = 32;
        int gpu_local_max = sqrt(opencl_context.max_thread_block_size());
        if (gpu_local_max < thread_block_size) {
          thread_block_size = gpu_local_max;
        }
        int parts = (inv.rows()+thread_block_size-1)/thread_block_size;
        try {
          opencl_kernels::lower_tri_inverse_step1(
              cl::NDRange(parts*thread_block_size),
              cl::NDRange(thread_block_size),
              inv.buffer(), inv.rows());
        } catch (cl::Error& e) {
          check_opencl_error("inverse step1", e);
        }
        parts /= 2;
        int result_matrix_dim = thread_block_size;
        while (parts > 0) {
          std::cout << "b: " << result_matrix_dim << std::endl;
          opencl_kernels::lower_tri_inverse_step2(
              cl::NDRange(result_matrix_dim, result_matrix_dim, parts),
              cl::NDRange(thread_block_size, thread_block_size, 1),
              inv.buffer(), temp.buffer(), inv.rows(), result_matrix_dim, result_matrix_dim );
          opencl_kernels::lower_tri_inverse_step3(
              cl::NDRange(result_matrix_dim, result_matrix_dim, parts),
              cl::NDRange(thread_block_size, thread_block_size, 1),
              inv.buffer(), temp.buffer(), inv.rows(), result_matrix_dim, result_matrix_dim );
          parts /= 2;
          result_matrix_dim *= 2;        
        }
      } catch (const cl::Error& e) {
        check_opencl_error("inverse", e);
      }
      return inv;
    }
  }
}

#endif
#endif
