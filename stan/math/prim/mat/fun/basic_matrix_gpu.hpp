#ifndef STAN_MATH_PRIM_ARR_BASIC_MATRIX_GPU_HPP
#define STAN_MATH_PRIM_ARR_BASIC_MATRIX_GPU_HPP

#include <stan/math/prim/mat/fun/ocl_gpu.hpp>
#include <stan/math/prim/arr/fun/matrix_gpu.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <iostream>
#include <string>
#include <map>

/** @file stan/math/prim/mat/fun/basic_matrix_gpu.hpp
*    @brief basic_matrix - basic matrix operations:
*    matrix transpose, copy matrices on the GPU,  ...
*/

namespace stan {
  namespace math {

    //transpose the matrix B and store it in matrix A
    void transpose(stan::math::matrix_gpu & dst, stan::math::matrix_gpu & src) {
      if(!(dst.rows == src.cols && dst.cols == src.rows)) {
        app_error("If the input matrix is MxN, the output matrix should be NxM.");
      }

      cl::Kernel kernel = stan::math::get_kernel("transpose");
      cl::CommandQueue cmdQueue = stan::math::get_queue();

      try {
        kernel.setArg(0, dst.buffer());
        kernel.setArg(1, src.buffer());
        kernel.setArg(2, src.rows);
        kernel.setArg(3, src.cols);
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(src.rows, src.cols),
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
      stan::math::matrix_gpu temp(A.cols, A.rows);
      stan::math::transpose(temp, A);
      int dim_temp = A.rows;
      A.rows = A.cols;
      A.cols = dim_temp;
      stan::math::copy(A, temp);

    }

    void zeros(stan::math::matrix_gpu & A, triangularity part = NONE) {
      cl::Kernel kernel = stan::math::get_kernel("zeros");
      cl::CommandQueue cmdQueue = stan::math::get_queue();
      try {
        kernel.setArg(0, A.buffer());
        kernel.setArg(1, A.rows);
        kernel.setArg(2, A.cols);
        kernel.setArg(3, part);
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(A.rows, A.cols),
          cl::NullRange,
          NULL,
          NULL);
      } catch (cl::Error& e) {
        check_ocl_error(e);
      }
    }

    void identity(stan::math::matrix_gpu & A) {
      if (A.rows!= A.cols) {
        app_error("Cannot make an identity matrix in a non-square matrix!");
      }
      cl::Kernel kernel = stan::math::get_kernel("identity");
      cl::CommandQueue cmdQueue = stan::math::get_queue();

      try {
        kernel.setArg(0, A.buffer());
        kernel.setArg(1, A.rows);
        kernel.setArg(2, A.cols);
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(A.rows, A.cols),
          cl::NullRange,
          NULL,
          NULL);
      } catch (cl::Error& e) {
        check_ocl_error(e);
      }
    }

    void copy_triangular(stan::math::matrix_gpu & src,
     stan::math::matrix_gpu & dst, triangularity lower_upper) {
      if (!(src.rows == dst.rows && src.cols == dst.cols)) {
        app_error("The dimensions of the input and output matrices should match.");
      }
      cl::Kernel kernel = stan::math::get_kernel("copy_triangular");
      cl::CommandQueue cmdQueue = stan::math::get_queue();
      try {
        kernel.setArg(0, src.buffer());
        kernel.setArg(1, dst.buffer());
        kernel.setArg(2, dst.rows);
        kernel.setArg(3, dst.cols);
        kernel.setArg(4, lower_upper);
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(dst.rows, dst.cols),
          cl::NullRange,
          NULL,
          NULL);
      } catch (cl::Error& e) {
        check_ocl_error(e);
      }
    }

    void copy_triangular_transposed(stan::math::matrix_gpu & A,
     copy_transposed_triangular lower_upper) {

      if (!(A.rows == A.cols)) {
        app_error("The matrix in the triangular"
        "transposed copy is non-square.");
      }
      cl::Kernel kernel = stan::math::get_kernel("copy_triangular_transposed");
      cl::CommandQueue cmdQueue = stan::math::get_queue();
      try {
        kernel.setArg(0, A.buffer());
        kernel.setArg(1, A.rows);
        kernel.setArg(2, A.cols);
        kernel.setArg(3, lower_upper);
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(A.rows, A.cols),
          cl::NullRange,
          NULL,
          NULL);
      } catch (cl::Error& e) {
        check_ocl_error(e);
      }
    }

    void add(stan::math::matrix_gpu & C,
     stan::math::matrix_gpu & A, stan::math::matrix_gpu & B) {

      if (!( A.rows == B.rows && C.rows == B.rows && A.cols == B.cols && C.cols == B.cols )) {
        app_error("The matrix dimensions in matrix addition do not match!");
      }
      cl::Kernel kernel = stan::math::get_kernel("add");
      cl::CommandQueue cmdQueue = stan::math::get_queue();
      try {
        kernel.setArg(0, C.buffer());
        kernel.setArg(1, A.buffer());
        kernel.setArg(2, B.buffer());
        kernel.setArg(3, A.rows);
        kernel.setArg(4, A.cols);
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(A.rows, A.cols),
          cl::NullRange,
          NULL,
          NULL);
      } catch (cl::Error& e) {
        check_ocl_error(e);
      }
    }

    void subtract(stan::math::matrix_gpu & C, stan::math::matrix_gpu & A,
     stan::math::matrix_gpu & B) {

      if (!( A.rows == B.rows && C.rows == B.rows && A.cols == B.cols && C.cols == B.cols)) {
        app_error("The matrix dimensions in matrix addition do not match!");
      }
      cl::Kernel kernel = stan::math::get_kernel("subtract");
      cl::CommandQueue cmdQueue = stan::math::get_queue();

      try {
        kernel.setArg(0, C.buffer());
        kernel.setArg(1, A.buffer());
        kernel.setArg(2, B.buffer());
        kernel.setArg(3, A.rows);
        kernel.setArg(4, A.cols);
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(A.rows, A.cols),
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
