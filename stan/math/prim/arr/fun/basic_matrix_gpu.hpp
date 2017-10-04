#ifndef STAN_MATH_PRIM_ARR_BASIC_MATRIX_GPU_HPP
#define STAN_MATH_PRIM_ARR_BASIC_MATRIX_GPU_HPP

#include <stan/math/prim/mat/fun/ocl.hpp>
#include <stan/math/prim/mat/fun/matrix.hpp>
#include <iostream>
#include <string>
#include <map>
#include <Eigen/Dense>

/** @file stan/math/prim/mat/fun/basic_matrix.hpp
*    @brief basic_matrix - basic matrix operations:
*    matrix transpose, copy matrices on the GPU,  ...
*/

namespace stan {
  namespace math {

    //transpose the matrix B and store it in matrix A
    void transpose(stan::math::matrix_gpu & dst, stan::math::matrix_gpu & src) {
      if(!(dst.M == src.N && dst.N == src.M)) {
        app_error("output matrix dimensions in matrix transpose!\n"
        "\nIf the input matrix is MxN, the output matrix should be NxM.");
      }

      cl::Kernel kernel = stan::math::get_kernel("transpose");
      cl::CommandQueue cmdQueue = stan::math::get_queue();

      try {
        kernel.setArg(0, dst.buffer());
        kernel.setArg(1, src.buffer());
        kernel.setArg(2, src.M);
        kernel.setArg(3, src.N);
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(src.M, src.N),
          cl::NullRange,
          NULL,
          NULL);
      } catch (cl::Error& e) {
        check_ocl_error(e);
      }
    }

    //when transpose with a single matrix is called,
    // a temporary matrix is created and used to transpose
    void transpose(stan::math::matrix_gpu & A) {
      stan::math::matrix_gpu temp(A.N, A.M);
      stan::math::transpose(temp, A);
      int dim_temp = A.M;
      A.M = A.N;
      A.N = dim_temp;
      stan::math::copy(A, temp);

    }

    void zeros(stan::math::matrix_gpu & A, triangularity part = NONE) {
      cl::Kernel kernel = stan::math::get_kernel("zeros");
      cl::CommandQueue cmdQueue = stan::math::get_queue();
      try {
        kernel.setArg(0, A.buffer());
        kernel.setArg(1, A.M);
        kernel.setArg(2, A.N);
        kernel.setArg(3, part);
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(A.M, A.N),
          cl::NullRange,
          NULL,
          NULL);
      } catch (cl::Error& e) {
        check_ocl_error(e);
      }
    }

    void identity(stan::math::matrix_gpu & A) {
      if(A.M! = A.N) {
        app_error("Can not make an identity matrix in a non-square matrix!");
      }
      cl::Kernel kernel = stan::math::get_kernel("identity");
      cl::CommandQueue cmdQueue = stan::math::get_queue();

      try {
        kernel.setArg(0, A.buffer());
        kernel.setArg(1, A.M);
        kernel.setArg(2, A.N);
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(A.M, A.N),
          cl::NullRange,
          NULL,
          NULL);
      } catch (cl::Error& e) {
        check_ocl_error(e);
      }
    }

    void copy_triangular(stan::math::matrix_gpu & src,
     stan::math::matrix_gpu & dst, triangularity lower_upper) {
      if(!(src.M == dst.M && src.N == dst.N)) {
        app_error("output matrix dimensions in lower triangular copy error!\n"
        "\nThe dimensions of the input and output matrices should match.");
      }
      cl::Kernel kernel = stan::math::get_kernel("copy_triangular");
      cl::CommandQueue cmdQueue = stan::math::get_queue();
      try {
        kernel.setArg(0, src.buffer());
        kernel.setArg(1, dst.buffer());
        kernel.setArg(2, dst.M);
        kernel.setArg(3, dst.N);
        kernel.setArg(4, lower_upper);
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(dst.M, dst.N),
          cl::NullRange,
          NULL,
          NULL);
      } catch (cl::Error& e) {
        check_ocl_error(e);
      }
    }

    void copy_triangular_transposed(stan::math::matrix_gpu & A,
     copy_transposed_triangular lower_upper) {

      if(!(A.M == A.N)) {
        app_error("The matrix in the triangular"
        "transposed copy is non-square.");
      }
      cl::Kernel kernel = stan::math::get_kernel("copy_triangular_transposed");
      cl::CommandQueue cmdQueue = stan::math::get_queue();
      try {
        kernel.setArg(0, A.buffer());
        kernel.setArg(1, A.M);
        kernel.setArg(2, A.N);
        kernel.setArg(3, lower_upper);
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(A.M, A.N),
          cl::NullRange,
          NULL,
          NULL);
      } catch (cl::Error& e) {
        check_ocl_error(e);
      }
    }

    void add(stan::math::matrix_gpu & C,
     stan::math::matrix_gpu & A, stan::math::matrix_gpu & B) {

      if(!( A.M == B.M && C.M == B.M && A.N == B.N && C.N == B.N )) {
        app_error("The matrix dimensions in matrix addition do not match!");
      }
      cl::Kernel kernel = stan::math::get_kernel("add");
      cl::CommandQueue cmdQueue = stan::math::get_queue();
      try {
        kernel.setArg(0, C.buffer());
        kernel.setArg(1, A.buffer());
        kernel.setArg(2, B.buffer());
        kernel.setArg(3, A.M);
        kernel.setArg(4, A.N);
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(A.M, A.N),
          cl::NullRange,
          NULL,
          NULL);
      } catch (cl::Error& e) {
        check_ocl_error(e);
      }
    }

    void subtract(stan::math::matrix_gpu & C, stan::math::matrix_gpu & A,
     stan::math::matrix_gpu & B) {

      if(!( A.M == B.M && C.M == B.M && A.N == B.N && C.N == B.N )) {
        app_error("The matrix dimensions in matrix addition do not match!");
      }
      cl::Kernel kernel = stan::math::get_kernel("subtract");
      cl::CommandQueue cmdQueue = stan::math::get_queue();

      try {
        kernel.setArg(0, C.buffer());
        kernel.setArg(1, A.buffer());
        kernel.setArg(2, B.buffer());
        kernel.setArg(3, A.M);
        kernel.setArg(4, A.N);
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(A.M, A.N),
          cl::NullRange,
          NULL,
          NULL);
      } catch (cl::Error& e) {
        check_ocl_error(e);
      }
    }
  }
}

#endif
