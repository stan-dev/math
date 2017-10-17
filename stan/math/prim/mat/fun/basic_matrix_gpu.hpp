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
    void transpose(matrix_gpu & dst, matrix_gpu & src) {
      if ((dst.rows() != src.cols() || dst.cols() != src.rows())) {
        app_error("If the input matrix is MxN, the output matrix should be NxM.");
      }

      cl::Kernel kernel = get_kernel("transpose");
      cl::CommandQueue cmdQueue = get_queue();

      try {
        kernel.setArg(0, dst.buffer());
        kernel.setArg(1, src.buffer());
        kernel.setArg(2, src.rows());
        kernel.setArg(3, src.cols());
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(src.rows(), src.cols()),
          cl::NullRange,
          NULL,
          NULL);
      } catch (const cl::Error& e) {
        check_ocl_error(e);
      }
    }

    //when transpose with a single matrix is called,
    // a temporary matrix is created and used to transpose
    void transpose(matrix_gpu & A) {
      matrix_gpu temp(A.cols(), A.rows());
      transpose(temp, A);
      int dim_temp = A.rows();
      A.rows_ = A.cols();
      A.cols_ = dim_temp;
      copy(A, temp);

    }

    void zeros(matrix_gpu & A, triangularity part = NONE) {
      cl::Kernel kernel = get_kernel("zeros");
      cl::CommandQueue cmdQueue = get_queue();
      try {
        kernel.setArg(0, A.buffer());
        kernel.setArg(1, A.rows());
        kernel.setArg(2, A.cols());
        kernel.setArg(3, part);
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(A.rows(), A.cols()),
          cl::NullRange,
          NULL,
          NULL);
      } catch (const cl::Error& e) {
        check_ocl_error(e);
      }
    }

    void identity(matrix_gpu & A) {
      if (A.rows()!= A.cols()) {
        app_error("Cannot make an identity matrix in a non-square matrix!");
      }
      cl::Kernel kernel = get_kernel("identity");
      cl::CommandQueue cmdQueue = get_queue();

      try {
        kernel.setArg(0, A.buffer());
        kernel.setArg(1, A.rows());
        kernel.setArg(2, A.cols());
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(A.rows(), A.cols()),
          cl::NullRange,
          NULL,
          NULL);
      } catch (const cl::Error& e) {
        check_ocl_error(e);
      }
    }

    void copy_triangular(matrix_gpu & src,
     matrix_gpu & dst, triangularity lower_upper) {
      if ((src.rows() != dst.rows() || src.cols() != dst.cols())) {
        app_error("The dimensions of the input and output matrices should match.");
      }
      cl::Kernel kernel = get_kernel("copy_triangular");
      cl::CommandQueue cmdQueue = get_queue();
      try {
        kernel.setArg(0, src.buffer());
        kernel.setArg(1, dst.buffer());
        kernel.setArg(2, dst.rows());
        kernel.setArg(3, dst.cols());
        kernel.setArg(4, lower_upper);
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(dst.rows(), dst.cols()),
          cl::NullRange,
          NULL,
          NULL);
      } catch (const cl::Error& e) {
        check_ocl_error(e);
      }
    }

    void copy_triangular_transposed(matrix_gpu & A,
     copy_transposed_triangular lower_upper) {

      if (A.rows() != A.cols()) {
        app_error("The matrix in the triangular"
        "transposed copy is non-square.");
      }
      cl::Kernel kernel = get_kernel("copy_triangular_transposed");
      cl::CommandQueue cmdQueue = get_queue();
      try {
        kernel.setArg(0, A.buffer());
        kernel.setArg(1, A.rows());
        kernel.setArg(2, A.cols());
        kernel.setArg(3, lower_upper);
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(A.rows(), A.cols()),
          cl::NullRange,
          NULL,
          NULL);
      } catch (const cl::Error& e) {
        check_ocl_error(e);
      }
    }

    void add(matrix_gpu & C,
     matrix_gpu & A, matrix_gpu & B) {

      if ( A.rows() != B.rows() || C.rows() != B.rows() || A.cols() != B.cols() || C.cols() != B.cols() ) {
        app_error("The matrix dimensions in matrix addition do not match!");
      }
      cl::Kernel kernel = get_kernel("add");
      cl::CommandQueue cmdQueue = get_queue();
      try {
        kernel.setArg(0, C.buffer());
        kernel.setArg(1, A.buffer());
        kernel.setArg(2, B.buffer());
        kernel.setArg(3, A.rows());
        kernel.setArg(4, A.cols());
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(A.rows(), A.cols()),
          cl::NullRange,
          NULL,
          NULL);
      } catch (const cl::Error& e) {
        check_ocl_error(e);
      }
    }

    void subtract(matrix_gpu & C, matrix_gpu & A,
     matrix_gpu & B) {

      if ( A.rows() != B.rows() || C.rows() != B.rows() || A.cols() != B.cols() || C.cols() != B.cols() ) {
        app_error("The matrix dimensions in matrix addition do not match!");
      }
      cl::Kernel kernel = get_kernel("subtract");
      cl::CommandQueue cmdQueue = get_queue();

      try {
        kernel.setArg(0, C.buffer());
        kernel.setArg(1, A.buffer());
        kernel.setArg(2, B.buffer());
        kernel.setArg(3, A.rows());
        kernel.setArg(4, A.cols());
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(A.rows(), A.cols()),
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
