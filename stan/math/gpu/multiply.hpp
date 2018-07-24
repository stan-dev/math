#ifndef STAN_MATH_GPU_MULTIPLY_HPP
#define STAN_MATH_GPU_MULTIPLY_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/matrix_gpu.hpp>
#include <Eigen/Dense>

namespace stan {
  namespace math {
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
        opencl_context.set_kernel_args(kernel, temp.buffer(),
                                   A.buffer(), scalar, A.rows(),
                                   A.cols());
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
     */
    inline matrix_gpu multiply(matrix_gpu & A, matrix_gpu & B) {
      check_size_match("multiply (GPU)", "A.cols()", A.cols(),
       "B.rows()", B.rows());
      matrix_gpu temp(A.rows(), B.cols());
      if ( temp.size() == 0 )
        return temp;
      cl::Kernel kernel = opencl_context.get_kernel("matrix_multiply");
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
        opencl_context.set_kernel_args(kernel, A.rows(), B.cols(),
                                   B.rows(), A.buffer(), B.buffer(),
                                   temp.buffer());
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
    }
  }
}

#endif
#endif
