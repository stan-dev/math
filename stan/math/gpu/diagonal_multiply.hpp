#ifndef STAN_MATH_GPU_DIAGONAL_MULTIPLY_HPP
#define STAN_MATH_GPU_DIAGONAL_MULTIPLY_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/matrix_gpu.hpp>
#include <Eigen/Dense>

namespace stan {
namespace math {
  /**
   * Multiplies the diagonal of the specified matrix on the GPU 
   * with the specified scalar.
   * 
   * @param A input matrix
   * @param scalar scalar
   * 
   * @return copy of the input matrix with the diagonal multiplied by scalar
   *      
   */
inline matrix_gpu diagonal_multiply(const matrix_gpu & A, const double scalar) {
    matrix_gpu B(A);
    if (B.size() == 0)
      return B;
    cl::Kernel kernel = opencl_context.get_kernel("scalar_mul_diagonal");
    cl::CommandQueue cmdQueue = opencl_context.queue();
    int min_dim = B.rows();
    if (B.cols() < min_dim)
      min_dim = B.cols();
    try {
      opencl_context.set_kernel_args(kernel, B.buffer(), scalar,
                                  B.rows(), B.cols());
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
    return B;
  }
}
}

#endif
#endif
