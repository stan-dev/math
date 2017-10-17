#ifndef STAN_MATH_PRIM_MAT_FUN_MATRIX_INVERSE_GPU_HPP
#define STAN_MATH_PRIM_MAT_FUN_MATRIX_INVERSE_GPU_HPP

#include <stan/math/prim/mat/fun/ocl_gpu.hpp>
#include <stan/math/prim/arr/fun/matrix_gpu.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <map>

/** @file stan/math/prim/mat/fun/inverse_gpu.hpp
*    @brief matrix_inverse - functions for matrix inversion:
*     lower triangular,  upper triangular,  regular,  ...
*/

//CURRENTLY ONLY SUPPORTS LOWER TRIANGULAR
namespace stan {
  namespace math {
    //transpose the matrix B and store it in matrix A
    void lower_triangular_inverse(matrix_gpu & A) {
      if (A.rows() != A.cols()) {
        app_error("the input matrix of the inverse is not square");
      }

      cl::Kernel kernel_step1 = get_kernel("lower_tri_inv_step1");
      cl::Kernel kernel_step2 = get_kernel("lower_tri_inv_step2");
      cl::Kernel kernel_step3 = get_kernel("lower_tri_inv_step3");
      cl::CommandQueue cmdQueue = get_queue();

      try{

        int parts = 32;
        //this will be managed by the library core with self-adaptive strategies
        if (A.rows() < 65)
          parts = 1;

        if (A.rows() > 2500)
          parts = 64;

        matrix_gpu temp(A.rows(),  A.cols() * 2);

        int remainder = A.rows() % parts;
        int part_size_fixed = A.rows()/parts;

        std::vector<int> stl_sizes (parts,part_size_fixed);
        
        for(int i = 0; i < remainder; i++) {
            stl_sizes[i]++;
        }
        cl::Context ctx = get_context();
        cl::Buffer sizes = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(int) * parts);

        cmdQueue.enqueueWriteBuffer(sizes, CL_TRUE, 0,
         sizeof(int) * parts, &stl_sizes[0]);

        kernel_step1.setArg(0, A.buffer());
        kernel_step1.setArg(1, temp.buffer());
        kernel_step1.setArg(2, remainder);
        kernel_step1.setArg(3, part_size_fixed);
        kernel_step1.setArg(4, A.rows());


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
          kernel_step2.setArg(6, A.rows());

          kernel_step3.setArg(0, A.buffer());
          kernel_step3.setArg(1, sizes);
          kernel_step3.setArg(2, temp.buffer());
          kernel_step3.setArg(3, repeat);
          kernel_step3.setArg(4, remainder);
          kernel_step3.setArg(5, part_size_fixed);
          kernel_step3.setArg(6, A.rows());

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
    }
  }
}

#endif
