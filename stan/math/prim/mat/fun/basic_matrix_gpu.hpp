#ifndef STAN_MATH_PRIM_ARR_BASIC_MATRIX_GPU_HPP
#define STAN_MATH_PRIM_ARR_BASIC_MATRIX_GPU_HPP

#include <stan/math/prim/mat/fun/ocl_gpu.hpp>
#include <stan/math/prim/arr/fun/matrix_gpu.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/mat/err/check_gpu.hpp>
#include <stan/math/prim/mat/err/check_matching_dims.hpp>
#include <iostream>
#include <string>
#include <map>

/** @file stan/math/prim/mat/fun/basic_matrix_gpu.hpp
*    @brief basic_matrix - basic matrix operations:
*    matrix transpose, copy matrices on the GPU,  ...
*/

namespace stan {
  namespace math {

    /**
     * Stores the tranpose of the second matrix
     * to the first matrix. Both matrices 
     * are on the  GPU.
     * 
     * @param dst the result matrix
     * @param src the input matrix
     * 
     * @throw <code>std::invalid_argument</code> if the 
     * dimensions of the matrices do not match     
     */
    void transpose(matrix_gpu & dst, matrix_gpu & src) {
      check_size_match("transpose (GPU)",
        "dst.rows()", dst.rows(), "src.cols()", src.cols());
      check_size_match("transpose (GPU)",
        "dst.cols()", dst.cols(), "src.rows()", src.rows());

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

    /**
     * Transposes the matrix on the GPU
     * 
     * @param A matrix
     * 
     */
    void transpose(matrix_gpu & A) {
      matrix_gpu temp(A.cols(), A.rows());
      transpose(temp, A);
      int dim_temp = A.rows();
      A.rows_ = A.cols();
      A.cols_ = dim_temp;
      stan::math::copy(A, temp); //NOLINT
    }

    /**
     * Stores zeros in the matrix on the GPU. 
     * It support writing the zeros to lower
     * triangular, upper triangular or the 
     * whole matrix.
     * 
     * @param A matrix
     * @param part optional parameter
     * that describes where to assign zeros: 
     * LOWER - lower triangular
     * UPPER - upper triangular
     * if the parameter is not specified, 
     * zeros are assigned to the whole matrix.
     * 
     */
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

    /**
     * Stores the identity matrix to the
     * assigned matrix on the GPU
     * 
     * @param A matrix
     * 
     * @throw <code>std::invalid_argument</code> if the 
     * matrix is not square 
     * 
     */
    void identity(matrix_gpu & A) {
      check_square("identity (GPU)", "A", A);
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

    /**
     * Copies the lower or upper 
     * triangular of the source matrix to the 
     * destination matrix. Both matrices are
     * stored on the GPU. 
     * 
     * @param src the source matrix
     * @param dst the destination matrix
     * @param lower_upper enum to describe
     * which part of the matrix to copy: 
     * LOWER - copies the lower triangular
     * UPPER - copes the upper triangular
     * 
     * @throw <code>std::invalid_argument</code> if the 
     * matrices do not have matching dimensions
     * 
     */
    void copy_triangular(matrix_gpu & src,
     matrix_gpu & dst, triangularity lower_upper) {
      check_matching_dims("copy_triangular (GPU)", "src", src, "dst", dst);
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

    /**
     * Copies the lower triangular 
     * to the upper triangular of the 
     * matrix on the GPU.
     * 
     * @param A matrix
     * @param lower_upper enum to describe
     * which copy operation to perform 
     * LOWER_TO_UPPER_TRIANGULAR - copies the lower 
     *   triangular to the upper triangular
     * UPPER_TO_LOWER_TRIANGULAR - copes the upper
     *  triangular to the lower triangular
     * 
     * @throw <code>std::invalid_argument</code> if the 
     * matrices do not have matching dimensions
     * 
     */
    void copy_triangular_transposed(matrix_gpu & A,
     copy_transposed_triangular lower_upper) {
      check_square("copy_triangular_transposed (GPU)", "A", A);
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
    /**
     * Matrix addition on the GPU
     * Adds the first two matrices and stores 
     * the result in the third matrix (C=A+B)
     * 
     * @param A first matrix
     * @param B second matrix
     * @param C result matrix
     * 
     * @throw <code>std::invalid_argument</code> if the 
     * matrices do not have matching dimensions
     * 
     */
    void add(matrix_gpu& A,
     matrix_gpu& B, matrix_gpu& C) {
      check_matching_dims("add (GPU)", "A", A, "B", B);
      check_matching_dims("add (GPU)", "B", B, "C", C);
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

    /**
     * Matrix subtraction on the GPU
     * Subtracts the second matrix
     * from the first matrix and stores 
     * the result in the third matrix (C=A-B)
     * 
     * @param A first matrix
     * @param B second matrix
     * @param C result matrix
     * 
     * @throw <code>std::invalid_argument</code> if the 
     * matrices do not have matching dimensions
     * 
     */
    void subtract(matrix_gpu & C, matrix_gpu & A,
     matrix_gpu & B) {
      check_matching_dims("subtract (GPU)", "A", A, "B", B);
      check_matching_dims("subtract (GPU)", "B", B, "C", C);
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

