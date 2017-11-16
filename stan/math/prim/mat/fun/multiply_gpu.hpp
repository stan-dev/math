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
    /**
     * Multiplies the diagonal of the specified matrix on the GPU 
     * with the specified scalar.
     * 
     * @param[in,out] A matrix
     * @param scalar scalar
     *      
     */
    void diagonal_multiply(matrix_gpu & A, double scalar) {
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
        check_ocl_error("diagonal_multiply", e);
      }
    }

    /**
     * Multiplies the specified matrix on the GPU 
     * with the specified scalar.
     * 
     * @param A matrix
     * @param scalar scalar
     * 
     * @return matrix multipled with scalar
     *      
     */
    matrix_gpu multiply(matrix_gpu & A,  double scalar) {
      matrix_gpu temp(A.rows(), A.cols());
      if (A.size() == 0)
        return temp;
      cl::Kernel kernel = get_kernel("scalar_mul");
      cl::CommandQueue cmdQueue = get_queue();
      try {
        kernel.setArg(0, temp.buffer());
        kernel.setArg(1, A.buffer());
        kernel.setArg(2, scalar);
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
        check_ocl_error("multiply scalar", e);
      }
      return temp;
    }

    /**
     * Multiplies the specified matrix on the GPU 
     * with the specified scalar.
     * 
     * @param scalar scalar
     * @param A matrix     
     * 
     * @return matrix multipled with scalar
     *      
     */
    matrix_gpu multiply(double scalar, matrix_gpu & A) {
      matrix_gpu temp(A.rows(), A.cols());
      if (A.size() == 0)
        return temp;
      cl::Kernel kernel = get_kernel("scalar_mul");
      cl::CommandQueue cmdQueue = get_queue();
      try {
        kernel.setArg(0, temp.buffer());
        kernel.setArg(1, A.buffer());
        kernel.setArg(2, scalar);
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
        check_ocl_error("multiply scalar", e);
      }
      return temp;
    }

    /**
     * Computes the product of the first two specified GPU matrices.  The number of
     * columns in the first matrix must be the same as the number of rows
     * in the second matrix. The number of rows of the resulting matrix must be the
     * same as the number of rows in the first matrix and the number of columnds
     * of the resulting matrix must be the same as the number of columns of the 
     * second matrix. 
     * 
     * Computes the matrix multiplication C = A x B
     * 
     * @param A first matrix
     * @param B second matrix
     * 
     * @return the product of the first and second matrix
     * 
     * @throw <code>std::invalid_argument</code> if the 
     *   number of columns in A and rows in B do not match
     */
    matrix_gpu multiply(matrix_gpu & A, matrix_gpu & B) {
      check_size_match("multiply (GPU)", "A.cols()", A.cols(),
       "B.rows()", B.rows());
      matrix_gpu temp(A.rows(), B.cols());
      cl::Kernel kernel = get_kernel("basic_multiply");
      cl::CommandQueue& cmdQueue = get_queue();
      try {
        int local = 8;
        std::cout << A.rows << " x " << B.cols() << std::endl;
        int Mpad = ((A.rows() + local-1)/local)*local;
        int Npad = ((B.cols() + local-1)/local)*local;
        kernel.setArg(0, A.rows());
        kernel.setArg(1, B.cols());
        kernel.setArg(2, B.rows());
        kernel.setArg(3, A.buffer());
        kernel.setArg(4, B.buffer());
        kernel.setArg(5, temp.buffer());
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(Mpad,  Npad),
          cl::NDRange(local, local),
          NULL,
          NULL);
      } catch (cl::Error& e) {
        check_ocl_error("multiply", e);
      }
      return temp;
    }
  }
}

#endif
