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
      /*cl::Kernel kernel_step1 =
        opencl_context.get_kernel("lower_tri_inverse_step1");
      cl::Kernel kernel_step2 =
        opencl_context.get_kernel("lower_tri_inverse_step2");
      cl::Kernel kernel_step3 =
        opencl_context.get_kernel("lower_tri_inverse_step3");*/
      cl::CommandQueue cmdQueue = opencl_context.queue();
      int parts = 1;
      if (inv.rows() < 65)
        parts = 1;
      if (inv.rows() > 2500)
        parts = 64;
      try {
        matrix_gpu temp(inv.rows(),  inv.cols() * 2);

        int remainder = inv.rows() % parts;
        int part_size_fixed = inv.rows()/parts;
        
        std::vector<int> stl_sizes(parts, part_size_fixed);

        for (int i = 0; i < remainder; i++) {
            stl_sizes[i]++;
        }
        cl::Context& ctx = opencl_context.context();
        cl::Buffer sizes = cl::Buffer(ctx, CL_MEM_READ_WRITE,
         sizeof(int) * parts);

        cmdQueue.enqueueWriteBuffer(sizes, CL_TRUE, 0,
         sizeof(int) * parts, &stl_sizes[0]);
        try {
          opencl_kernels::lower_tri_inverse_step1(
              cl::NDRange(parts), cl::NDRange(1),
              inv.buffer(), temp.buffer(), remainder, part_size_fixed, inv.rows());
        } catch (cl::Error& e) {
          check_opencl_error("multiply self transpose", e);
        }

        int repeat = 1;
        int sizePad;

        int local = 32;
        int gpu_local_max = sqrt(opencl_context.max_thread_block_size());
        if (gpu_local_max < local) {
          local = gpu_local_max;
        }
        int wpt = 4;

        for (int pp = parts; pp > 1; pp /= 2) {
          std::cout << "abc" << std::endl;
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
            check_opencl_error("multiply self transpose", e);
          }
          repeat *= 2;
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
