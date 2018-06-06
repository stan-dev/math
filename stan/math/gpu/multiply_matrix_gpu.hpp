#ifndef STAN_MATH_PRIM_MAT_FUN_MULTIPY_GPU_HPP
#define STAN_MATH_PRIM_MAT_FUN_MULTIPY_GPU_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/basic_matrix_gpu.hpp>
#include <Eigen/Dense>

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
    inline void diagonal_multiply(matrix_gpu & A, double scalar) {
      cl::Kernel kernel = opencl_context.get_kernel("scalar_mul_diagonal");
      cl::CommandQueue cmdQueue = opencl_context.queue();
      int min_dim = A.rows();
      if (A.cols() < min_dim)
        min_dim = A.cols();
      try {
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
        check_opencl_error("diagonal_multiply", e);
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
    inline matrix_gpu multiply(matrix_gpu & A,  double scalar) {
      matrix_gpu temp(A.rows(), A.cols());
      if (A.size() == 0)
        return temp;
      cl::Kernel kernel = opencl_context.get_kernel("scalar_mul");
      cl::CommandQueue cmdQueue = opencl_context.queue();
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
        check_opencl_error("multiply scalar", e);
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
    inline matrix_gpu multiply(double scalar, matrix_gpu & A) {
      return multiply(A, scalar);
    }

    /**
     * Computes the product of the two specified GPU matrices.  The number of
     * columns in the first matrix must be the same as the number of rows
     * in the second matrix. The number of rows of the resulting matrix must be the
     * same as the number of rows in the first matrix and the number of columns
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
     *//*
    inline matrix_gpu multiply(matrix_gpu & A, matrix_gpu & B) {
      check_size_match("multiply (GPU)", "A.cols()", A.cols(),
       "B.rows()", B.rows());
      matrix_gpu temp(A.rows(), B.cols());
      if ( temp.size() == 0 )
        return temp;
      cl::Kernel kernel = opencl_context.get_kernel("basic_multiply");
      cl::CommandQueue& cmdQueue = opencl_context.queue();
      int local = 32;
      int gpu_local_max = sqrt(opencl_context.max_workgroup_size());
      if (gpu_local_max < local) {
        local = gpu_local_max;
      }
      int wpt = 4;
      int Mpad = ((A.rows() + local-1)/local)*local;
      int Npad = ((B.cols() + local-1)/local)*local;
      try {
        kernel.setArg(0, A.rows());
        kernel.setArg(1, B.cols());
        kernel.setArg(2, B.rows());
        kernel.setArg(3, A.buffer());
        kernel.setArg(4, B.buffer());
        kernel.setArg(5, temp.buffer());
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(Mpad,  Npad/wpt),
          cl::NDRange(local, local/wpt),
          NULL,
          NULL);
      } catch (cl::Error& e) {
        check_opencl_error("multiply", e);
      }
      return temp;
    }*/
    /**
     * Computes the product of the specified GPU matrix with its transpose. 
	 * The input matrix must me squared and the output matrix must be the
	 * same size as the input matrix.
     * 
     * Computes the matrix multiplication C = A x A^T
     * 
     * @param A matrix
     * 
     * @return the product of the matrix and its transpose
     * 
     *//*
    inline matrix_gpu multiply_with_self_transposed(matrix_gpu & A) {
      // TODO(Rok): added a check if A is square
      matrix_gpu temp(A.rows(), A.rows());
      if (temp.size() == 0)
        return temp;
      cl::Kernel kernel = opencl_context.get_kernel("multiply_self_transposed");
      cl::CommandQueue& cmdQueue = opencl_context.queue();
      matrix_gpu AT = stan::math::transpose(A);
      int local = 32;
      int gpu_local_max = sqrt(opencl_context.max_workgroup_size());
      if (gpu_local_max < local) {
        local = gpu_local_max;
      }
      int wpt = 4;
      int Mpad = ((A.rows() + local-1)/local)*local;
      int Npad = ((AT.cols() + local-1)/local)*local;
      try {
        kernel.setArg(0, A.rows());
        kernel.setArg(1, AT.cols());
        kernel.setArg(2, AT.rows());
        kernel.setArg(3, A.buffer());
        kernel.setArg(4, AT.buffer());
        kernel.setArg(5, temp.buffer());
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(Mpad,  Npad/wpt),
          cl::NDRange(local, local/wpt),
          NULL,
          NULL);
      } catch (cl::Error& e) {
        check_opencl_error("multiply with self transposed", e);
      }
      return temp;
    }*/
	/**
     * Computes the lower triangular part of the product of the two specified
     * GPU matrices, where the first matrix is lower triangular. The number of
     * columns in the first matrix must be the same as the number of rows
     * in the second matrix. The number of rows of the resulting matrix must be the
     * same as the number of rows in the first matrix and the number of columns
     * of the resulting matrix must be the same as the number of columns of the 
     * second matrix. 
     * 
     * Computes the matrix multiplication L_out = L x B
     * 
     * @param L first matrix, which is lower triangular
     * @param B second matrix
     * 
     * @return lower triangular part of the product of the first and second matrix
     * 
     * @throw <code>std::invalid_argument</code> if the 
     *   number of columns in L and rows in B do not match
     *//*
    inline matrix_gpu multiply_lower_triangular(matrix_gpu & L,
        matrix_gpu & B) {
      check_size_match("lower triangular multiply (GPU)",
          "L.cols()", L.cols(), "B.rows()", B.rows());
      matrix_gpu temp(L.rows(), L.rows());
      if (temp.size() == 0)
        return temp;
      cl::Kernel kernel =
        opencl_context.get_kernel("multiply_lower_triangular");
      cl::CommandQueue& cmdQueue = opencl_context.queue();
      int local = 32;
      int gpu_local_max = sqrt(opencl_context.max_workgroup_size());
      if (gpu_local_max < local) {
        local = gpu_local_max;
      }
      int wpt = 4;
      int Mpad = ((L.rows() + local-1)/local)*local;
      int Npad = ((B.cols() + local-1)/local)*local;
      try {
        kernel.setArg(0, L.rows());
        kernel.setArg(1, B.cols());
        kernel.setArg(2, B.rows());
        kernel.setArg(3, L.buffer());
        kernel.setArg(4, B.buffer());
        kernel.setArg(5, temp.buffer());
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(Mpad,  Npad/wpt),
          cl::NDRange(local, local/wpt),
          NULL,
          NULL);
      } catch (cl::Error& e) {
        check_opencl_error("lower triangular multiply", e);
      }
      return temp;
    }*/
  }
}

#endif
#endif
