#ifndef STAN_MATH_GPU_SUBTRACT_OPENCL_HPP
#define STAN_MATH_GPU_SUBTRACT_OPENCL_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/matrix_gpu.hpp>
#include <stan/math/gpu/err/check_matching_dims.hpp>
#include <CL/cl.hpp>

/** @file stan/math/gpu/basic_matrix_gpu.hpp
*    @brief basic_matrix_gpu - basic matrix operations:
*    copy, copy lower/upper triangular, copy triangular transposed,
*    copy submatrix, init matrix with zeros, create an identity matrix,
*    add, subtract, transpose
*/

namespace stan {
  namespace math {

    /**
     * Matrix subtraction on the GPU
     * Subtracts the second matrix
     * from the first matrix and stores
     * the result in the third matrix (C=A-B)
     *
     * @param A first matrix
     * @param B second matrix
     *
     * @return subtraction result matrix
     *
     * @throw <code>std::invalid_argument</code> if the
     * input matrices do not have matching dimensions
     *
     */
    inline matrix_gpu subtract(matrix_gpu & A, matrix_gpu & B) {
      check_matching_dims("subtract (GPU)", "A", A, "B", B);
      matrix_gpu C(A.rows(), A.cols());
      if (A.size() == 0) {
        return C;
      }
      cl::Kernel kernel = opencl_context.get_kernel("subtract");
      cl::CommandQueue cmdQueue = opencl_context.queue();

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
        check_opencl_error("subtract", e);
      }
      return C;
    }
  }
}

#endif
#endif
