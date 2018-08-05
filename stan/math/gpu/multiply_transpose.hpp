#ifndef STAN_MATH_GPU_MULTIPLY_TRANSPOSE_HPP
#define STAN_MATH_GPU_MULTIPLY_TRANSPOSE_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/matrix_gpu.hpp>
#include <stan/math/gpu/err/check_square.hpp>
#include <Eigen/Dense>

namespace stan {
  namespace math {
    /**
     * Computes the product of a square GPU matrix with its transpose.
     *
     * Computes the matrix multiplication C = A x A^T
     *
     * @param A input matrix
     * @return the product of the input matrix and its transpose
     *
     */
    inline matrix_gpu multiply_transpose(const matrix_gpu & A) {
      matrix_gpu temp(A.rows(), A.rows());
      if (temp.size() == 0)
        return temp;
      // padding the matrices so the dimensions are divisible with local
      // improves performance becasuse we can omit if statements in the
      // multiply kernel
      int local = gpu::multiply_workgroup_size;
      int Mpad = ((A.rows() + local-1)/local)*local;
      int Npad = ((A.cols() + local-1)/local)*local;
      matrix_gpu tempPad(Mpad, Mpad);
      matrix_gpu Apad(Mpad, Npad);
      Apad.sub_block(A, 0, 0, 0, 0, A.rows(), A.cols());
      cl::Kernel kernel = opencl_context.get_kernel("multiply_transpose");
      cl::CommandQueue& cmdQueue = opencl_context.queue();
      int wpt = 4;
      try {
        opencl_context.set_kernel_args(kernel, Apad.buffer(), tempPad.buffer(),
                                   Apad.rows(), Apad.cols());
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
