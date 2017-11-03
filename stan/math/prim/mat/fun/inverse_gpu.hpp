#ifndef STAN_MATH_PRIM_MAT_FUN_MATRIX_INVERSE_GPU_HPP
#define STAN_MATH_PRIM_MAT_FUN_MATRIX_INVERSE_GPU_HPP

#include <stan/math/prim/mat/fun/ocl_gpu.hpp>
#include <stan/math/prim/arr/fun/matrix_gpu.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_gpu.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <map>
#include <vector>

/** @file stan/math/prim/mat/fun/inverse_gpu.hpp
*    @brief matrix_inverse - functions for matrix inversion:
*     lower triangular,  upper triangular,  regular,  ...
*/

// CURRENTLY ONLY SUPPORTS LOWER TRIANGULAR
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
    matrix_gpu lower_triangular_inverse(matrix_gpu & A) {
      check_square("lower_triangular_inverse (GPU)", "A", A);
      matrix_gpu inv(A.rows
      (), A.cols());
      copy(A, inv);//NOLINT
      cl::Kernel kernel_step1 = get_kernel("lower_tri_inv_step1");
      cl::Kernel kernel_step2 = get_kernel("lower_tri_inv_step2");
      cl::Kernel kernel_step3 = get_kernel("lower_tri_inv_step3");
      cl::CommandQueue cmdQueue = get_queue();

      try {
        int parts = 32;
        if (inv.rows() < 65)
          parts = 1;

        if (inv.rows() > 2500)
          parts = 64;

        matrix_gpu temp(inv.rows(),  inv.cols() * 2);

        int remainder = inv.rows() % parts;
        int part_size_fixed = inv.rows()/parts;

        std::vector<int> stl_sizes(parts, part_size_fixed);

        for (int i = 0; i < remainder; i++) {
            stl_sizes[i]++;
        }
        cl::Context ctx = get_context();
        cl::Buffer sizes = cl::Buffer(ctx, CL_MEM_READ_WRITE,
         sizeof(int) * parts);

        cmdQueue.enqueueWriteBuffer(sizes, CL_TRUE, 0,
         sizeof(int) * parts, &stl_sizes[0]);

        kernel_step1.setArg(0, inv.buffer());
        kernel_step1.setArg(1, temp.buffer());
        kernel_step1.setArg(2, remainder);
        kernel_step1.setArg(3, part_size_fixed);
        kernel_step1.setArg(4, inv.rows());


        cmdQueue.enqueueNDRangeKernel(
          kernel_step1,
          cl::NullRange,
          cl::NDRange(parts),
          cl::NDRange(1),
          NULL,
          NULL);

        int repeat = 1;
        int sizePad;
        for (int pp = parts; pp > 1; pp /= 2) {
          sizePad = (((part_size_fixed + 1) * repeat + 31) / 32) * 32;

          kernel_step2.setArg(0, inv.buffer());
          kernel_step2.setArg(1, sizes);
          kernel_step2.setArg(2, temp.buffer());
          kernel_step2.setArg(3, repeat);
          kernel_step2.setArg(4, remainder);
          kernel_step2.setArg(5, part_size_fixed);
          kernel_step2.setArg(6, inv.rows());

          kernel_step3.setArg(0, inv.buffer());
          kernel_step3.setArg(1, sizes);
          kernel_step3.setArg(2, temp.buffer());
          kernel_step3.setArg(3, repeat);
          kernel_step3.setArg(4, remainder);
          kernel_step3.setArg(5, part_size_fixed);
          kernel_step3.setArg(6, inv.rows());

           cmdQueue.enqueueNDRangeKernel(
          kernel_step2,
          cl::NullRange,
          cl::NDRange(sizePad, sizePad / 4, pp / 2),
          cl::NDRange(32, 8, 1),
          NULL,
          NULL);

            cmdQueue.enqueueNDRangeKernel(
          kernel_step3,
          cl::NullRange,
          cl::NDRange(sizePad, sizePad / 4, pp / 2),
          cl::NDRange(32, 8, 1),
          NULL,
          NULL);

          repeat *= 2;
        }
      } catch (const cl::Error& e) {
        check_ocl_error(e);
      }
      return inv;
    }
  }
}

#endif
