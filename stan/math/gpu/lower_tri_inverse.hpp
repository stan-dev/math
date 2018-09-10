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
      cl::CommandQueue cmdQueue = opencl_context.queue();
      int parts = (inv.rows()+1)/2;
      try {
        matrix_gpu temp((inv.rows()+1)/2,  (inv.cols()+1)/2);

        try {
          opencl_kernels::lower_tri_inverse_step1(
              cl::NDRange(parts),
              inv.buffer(), inv.rows());
        } catch (cl::Error& e) {
          check_opencl_error("inverse step1", e);
        }

        int local = 32;
        int gpu_local_max = sqrt(opencl_context.max_thread_block_size());
        if (gpu_local_max < local) {
          local = gpu_local_max;
        }

        parts /= 2;
        int block = 4;
        int result_matrix_dim = 2;
        std::cout << "max local: " << local << std::endl;
        while (parts > 0) {
        
          if (result_matrix_dim <= local){
            std::cout << result_matrix_dim << std::endl;
            std::cout << parts << std::endl;
            opencl_kernels::lower_tri_inverse_step2_small(
                cl::NDRange(result_matrix_dim, result_matrix_dim, parts),
                cl::NDRange(result_matrix_dim, result_matrix_dim, 1),
                inv.buffer(), result_matrix_dim, inv.rows());
          }
          parts /= 2;
          result_matrix_dim *= 2;
          block *= 2;
        
        }
        
        /*int wpt = 4;
        for (int pp = parts; pp > 1; pp /= 2) {
          sizePad = (((part_size_fixed + 1) * repeat + 31) / 32) * 32;
          try {
            opencl_kernels::lower_tri_inverse_step2(
                cl::NDRange(sizePad, sizePad / wpt, pp / 2),
                cl::NDRange(local, local/wpt, 1),
                inv.buffer(), sizes, temp.buffer(), repeat, remainder, part_size_fixed, inv.rows());
            opencl_kernels::lower_tri_inverse_step3(
                cl::NDRange(sizePad, sizePad / wpt, pp / 2),
                cl::NDRange(local, local/wpt, 1),
                inv.buffer(), sizes, temp.buffer(), repeat, remainder, part_size_fixed, inv.rows());
          } catch (cl::Error& e) {
            check_opencl_error("inverse step 2 and 3", e);
          }
          repeat *= 2;
        }*/
      } catch (const cl::Error& e) {
        check_opencl_error("inverse", e);
      }
      return inv;
    }
  }
}

#endif
#endif
