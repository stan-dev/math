#ifndef STAN_MATH_PRIM_MAT_FUN_MULTIPY_GPU_HPP
#define STAN_MATH_PRIM_MAT_FUN_MULTIPY_GPU_HPP

#include <stan/math/prim/mat/fun/ocl_gpu.hpp>
#include <stan/math/prim/arr/fun/matrix_gpu.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <map>

/** @file stan/math/prim/mat/fun/multiply_gpu.hpp
    @brief basic_matrix - basic matrix operations: matrix transpose,  copy matrices on the GPU,  ...
*/

namespace stan {
  namespace math {

    void diagonal_multiply_with_scalar(matrix_gpu & A,
     double scalar) {
      cl::Kernel kernel = get_kernel("scalar_mul_diagonal");
      cl::CommandQueue cmdQueue = get_queue();
      try {
        int min_dim = A.rows();
        if (A.cols() < min_dim)
          min_dim = A.cols();
        kernel.setArg(0, A.buffer());
        kernel.setArg(1, scalar),
        kernel.setArg(2, A.rows());
        kernel.setArg(3, A.cols());
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(min_dim),
          cl::NullRange,
          NULL,
          NULL);
      } catch (const cl::Error& e) {
        check_ocl_error(e);
      }
    }

    void multiply_with_scalar(matrix_gpu & A,  double scalar) {
      cl::Kernel kernel = get_kernel("scalar_mul");
      cl::CommandQueue cmdQueue = get_queue();
      try {
        kernel.setArg(0, A.buffer());
        kernel.setArg(1, scalar);
        kernel.setArg(2, A.rows());
        kernel.setArg(3, A.cols());
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

    void multiply(matrix_gpu & A, matrix_gpu & B,
     matrix_gpu & C ) {
      check_size_match("multiply (GPU)", "A.cols()", A.cols(),
       "B.rows()", B.rows());
      check_size_match("multiply (GPU)", "A.rows()", A.rows(),
       "C.rows()", C.rows());
      check_size_match("multiply (GPU)", "B.cols()", B.cols(),
       "C.cols()", C.cols());

      // Square matrices
      cl::Kernel kernel = get_kernel("basic_multiply");
      cl::CommandQueue cmdQueue = get_queue();
      try {
        int Mpad = ((A.rows() + 15)/16)*16;
        int Npad = ((B.cols() + 15)/16)*16;
        kernel.setArg(0, A.rows());
        kernel.setArg(1, B.cols());
        kernel.setArg(2, B.rows());
        kernel.setArg(3, A.buffer());
        kernel.setArg(4, B.buffer());
        kernel.setArg(5, C.buffer());
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(Mpad,  Npad),
          cl::NDRange(16, 16),
          NULL,
          NULL);
      } catch (cl::Error& e) {
        check_ocl_error(e);
      }
    }
  }
}

#endif
