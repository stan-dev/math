#ifndef STAN_MATH_GPU_MULTIPLY_SELF_TRANSPOSE_HPP
#define STAN_MATH_GPU_MULTIPLY_SELF_TRANSPOSE_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/matrix_gpu.hpp>
#include <stan/math/gpu/transpose.hpp>
#include <stan/math/gpu/err/check_square.hpp>
#include <Eigen/Dense>

namespace stan {
  namespace math {
    /**
     * Computes the product of the specified GPU matrix with its transpose. 
     * The input matrix must me squared and the output matrix must be the
	   * same size as the input matrix.
     * 
     * Computes the matrix multiplication C = A x A^T
     * 
     * @param A input matrix
     * 
     * @return the product of the input matrix and its transpose
     * 
     */
    inline matrix_gpu multiply_self_transpose(const matrix_gpu & A) {
      matrix_gpu temp(A.rows(), A.rows());
      if (temp.size() == 0)
        return temp;
      int local = 32;
      int gpu_local_max = sqrt(opencl_context.max_workgroup_size());
      if (gpu_local_max < local) {
        local = gpu_local_max;
      }
      int Mpad = ((A.rows() + local-1)/local)*local;
      int Npad = ((A.cols() + local-1)/local)*local;
      matrix_gpu tempPad(Mpad, Mpad);
      matrix_gpu Apad(Mpad, Npad);
      matrix_gpu ATpad(Npad, Mpad);
      matrix_gpu AT = stan::math::transpose(A);
      Apad.sub_block(A, 0, 0, 0, 0, A.rows(), A.cols());
      ATpad.sub_block(AT, 0, 0, 0, 0, AT.rows(), AT.cols());
      cl::Kernel kernel = opencl_context.get_kernel("multiply_self_transpose");
      cl::CommandQueue& cmdQueue = opencl_context.queue();
      int wpt = 4;
      try {
        opencl_context.set_kernel_args(kernel, Apad.rows(), ATpad.cols(),
                                   ATpad.rows(), Apad.buffer(), ATpad.buffer(),
                                   tempPad.buffer());
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(Mpad,  Mpad/wpt),
          cl::NDRange(local, local/wpt),
          NULL,
          NULL);
      } catch (cl::Error& e) {
        check_opencl_error("multiply self transpose", e);
      }
      temp.sub_block(tempPad, 0, 0, 0, 0, temp.rows(), temp.cols());      
      return temp;
    }
  }
}

#endif
#endif
