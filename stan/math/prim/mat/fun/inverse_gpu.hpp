#ifndef STAN_MATH_PRIM_MAT_FUN_MATRIX_INVERSE_GPU_HPP
#define STAN_MATH_PRIM_MAT_FUN_MATRIX_INVERSE_GPU_HPP

#include <stan/math/prim/mat/ocl_gpu.hpp>
#include <stan/math/prim/arr/fun/matrix_gpu.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <map>

/** @file stanmathcl/matrix_inverse.hpp
*    @brief matrix_inverse - functions for matrix inversion:
*     lower triangular,  upper triangular,  regular,  ...
*/

//CURRENTLY ONLY SUPPORTS LOWER TRIANGULAR
namespace stan {
  namespace math {
    //transpose the matrix B and store it in matrix A
    void lower_triangular_inverse(stanmathcl::Matrix & A) {
      if(!(A.M == A.N)) {
        app_error("the input matrix of the inverse is not square");
      }

      cl::Kernel kernel_step1 = stanmathcl::get_kernel("lower_tri_inv_step1");
      cl::Kernel kernel_step2 = stanmathcl::get_kernel("lower_tri_inv_step2");
      cl::Kernel kernel_step3 = stanmathcl::get_kernel("lower_tri_inv_step3");
      cl::CommandQueue cmdQueue = stanmathcl::get_queue();

      try{

        int parts = 32;
        //this will be managed by the library core with self-adaptive strategies
        if(A.M < 65)
          parts = 1;

        if(A.M > 2500)
          parts = 64;

        stanmathcl::Matrix temp(A.M,  A.N * 2);

        int remainder = A.M%parts;
        int part_size_fixed = A.M/parts;

        int* stl_sizes = new int[parts];

        for(int i = 0; i < parts; i++) {
          if(i < remainder)
            stl_sizes[i] = part_size_fixed+1;
          else
            stl_sizes[i] = part_size_fixed;
        }
        cl::Context ctx = stanmathcl::get_context();
        cl::Buffer sizes = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(int) * parts);

        cmdQueue.enqueueWriteBuffer(sizes, CL_TRUE, 0,
         sizeof(int) * parts, stl_sizes);

        kernel_step1.setArg(0, A.buffer());
        kernel_step1.setArg(1, temp.buffer());
        kernel_step1.setArg(2, remainder);
        kernel_step1.setArg(3, part_size_fixed);
        kernel_step1.setArg(4, A.M);


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

          kernel_step2.setArg(0, A.buffer());
          kernel_step2.setArg(1, sizes);
          kernel_step2.setArg(2, temp.buffer());
          kernel_step2.setArg(3, repeat);
          kernel_step2.setArg(4, remainder);
          kernel_step2.setArg(5, part_size_fixed);
          kernel_step2.setArg(6, A.M);

          kernel_step3.setArg(0, A.buffer());
          kernel_step3.setArg(1, sizes);
          kernel_step3.setArg(2, temp.buffer());
          kernel_step3.setArg(3, repeat);
          kernel_step3.setArg(4, remainder);
          kernel_step3.setArg(5, part_size_fixed);
          kernel_step3.setArg(6, A.M);

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

        delete stl_sizes;
      } catch (cl::Error& e) {

        check_ocl_error(e);

      }
    }
  }
}

#endif
