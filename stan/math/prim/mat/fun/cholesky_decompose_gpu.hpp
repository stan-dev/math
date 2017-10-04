#ifndef STAN_MATH_PRIM_MAT_FUN_CHOLESKY_DECOMPOSE_GPU_HPP
#define STAN_MATH_PRIM_MAT_FUN_CHOLESKY_DECOMPOSE_GPU_HPP

#include <stan/math/prim/mat/fun/ocl_gpu.hpp>
#include <stan/math/prim/arr/fun/matrix_gpu.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <map>

/*   @file stanmathcl/matrix_inverse.hpp
*    @brief matrix_inverse  -  functions for matrix inversion:
*     lower triangular,  upper triangular,  regular,  ...
*/

//CURRENTLY ONLY SUPPORTS LOWER TRIANGULAR
namespace stan {
  namespace math {

    template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    cholesky_gpu(stan::math::matrix_gpu & B) {

      stan::math::matrix_gpu A(B);
      if(!(A.M  ==  A.N)) {
        app_error("the input matrix of the decomposition is not square");
      }

      cl::Kernel kernel_chol_block = stanmathcl::get_kernel("cholesky_block");
      cl::Kernel kernel_left = stanmathcl::get_kernel("cholesky_left_update");
      cl::Kernel kernel_mid = stanmathcl::get_kernel("cholesky_mid_update");
      cl::Kernel kernel_zero = stanmathcl::get_kernel("cholesky_zero");
      cl::CommandQueue cmd_queue = stanmathcl::get_queue();

      try {
        cl::Context ctx = stanmathcl::get_context();
        //will be managed by the library core system
        int block = 64;
        int offset = 0;
        int local = 16;

        cl::Buffer buffer_V(ctx, CL_MEM_READ_WRITE,
         sizeof(T) * block * block * 4);
        cl::Buffer buffer_L(ctx, CL_MEM_READ_WRITE,
         sizeof(T) * block * A.M * 4);
        cl::Buffer buffer_D(ctx, CL_MEM_READ_WRITE,
         sizeof(T) * block * block * 4);

        kernel_chol_block.setArg(0, A.buffer());
        kernel_chol_block.setArg(1, offset);
        kernel_chol_block.setArg(2, A.M);
        kernel_chol_block.setArg(3, block);
        kernel_chol_block.setArg(4, buffer_V);
        kernel_chol_block.setArg(5, buffer_D);

        kernel_left.setArg(0, buffer_L);
        kernel_left.setArg(1, A.buffer());
        kernel_left.setArg(2, buffer_V);
        kernel_left.setArg(3, offset);
        kernel_left.setArg(4, block);
        kernel_left.setArg(5, A.M);

        kernel_mid.setArg(0, buffer_L);
        kernel_mid.setArg(1, A.buffer());
        kernel_mid.setArg(2, offset);
        kernel_mid.setArg(3, block);
        kernel_mid.setArg(4, A.M);

        kernel_zero.setArg(0, A.buffer());
        kernel_zero.setArg(1, A.M);

        int threadsLeft,  leftPad;
        int threadsMid,  midPad;
        while((offset + block) < (A.M - block)) {
          threadsLeft = A.M - offset - block;
          leftPad = ((threadsLeft + local - 1) / local) * local;
          threadsMid = A.M - offset - block;
          midPad = ((threadsMid + local - 1) / local) * local;
          kernel_left.setArg(3, offset);
          kernel_left.setArg(6, threadsLeft);
          kernel_mid.setArg(2, offset);
          kernel_mid.setArg(5, threadsMid);
          kernel_chol_block.setArg(1, offset);
          cmd_queue.enqueueNDRangeKernel(kernel_chol_block,
           cl::NullRange, cl::NDRange(block), cl::NDRange(block));
          cmd_queue.enqueueNDRangeKernel(kernel_left,
           cl::NullRange, cl::NDRange(leftPad, leftPad),
           cl::NDRange(local,  local));
          cmd_queue.enqueueNDRangeKernel(kernel_mid,
           cl::NullRange, cl::NDRange(midPad, midPad),
           cl::NDRange(local, local));
          offset += block;
        }

        int left = A.M - offset;
        if(left > 0) {
          kernel_chol_block.setArg(1, offset);
          kernel_chol_block.setArg(3, left);
          cmd_queue.enqueueNDRangeKernel(kernel_chol_block,
           cl::NullRange, cl::NDRange(left), cl::NDRange(left));
        }
      cmd_queue.enqueueNDRangeKernel(kernel_zero,
       cl::NullRange, cl::NDRange(A.M, A.M),  cl::NullRange);
      } catch (cl::Error& e) {
        check_ocl_error(e);
      }
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A_tmp;
      stan::math::copy(B, A_tmp);
      return A_tmp;
    }
  }
}

#endif
