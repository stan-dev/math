#ifndef STAN_MATH_GPU_COPY_TRIANGULAR_TRANSPOSED_OPENCL_HPP
#define STAN_MATH_GPU_COPY_TRIANGULAR_TRANSPOSED_OPENCL_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/matrix_gpu.hpp>
#include <stan/math/gpu/err/check_square.hpp>
#include <CL/cl.hpp>

namespace stan {
namespace math {
enum copy_transposed_triangular {
  LOWER_TO_UPPER_TRIANGULAR = 0,
  UPPER_TO_LOWER_TRIANGULAR = 1
};

/**
 * Copies a lower/upper triangular of a matrix to it's upper/lower.
 *
 * @param A matrix
 * @param lower_upper enum to describe
 * which copy operation to perform
 * LOWER_TO_UPPER_TRIANGULAR - copies the lower
 *   triangular to the upper triangular
 * UPPER_TO_LOWER_TRIANGULAR - copes the upper
 *  triangular to the lower triangular
 *
 * @throw <code>std::invalid_argument</code> if the matrix is not square.
 *
 */
inline void copy_triangular_transposed(matrix_gpu& A,
                                       copy_transposed_triangular lower_upper) {
  if (A.size() == 0) {
    return;
  }
  if (A.size() == 1) {
    return;
  }
  check_square("copy_triangular_transposed (GPU)", "A", A);
  cl::Kernel kernel = opencl_context.get_kernel("copy_triangular_transposed");
  cl::CommandQueue cmdQueue = opencl_context.queue();
  try {
    kernel.setArg(0, A.buffer());
    kernel.setArg(1, A.rows());
    kernel.setArg(2, A.cols());
    kernel.setArg(3, lower_upper);
    cmdQueue.enqueueNDRangeKernel(kernel, cl::NullRange,
                                  cl::NDRange(A.rows(), A.cols()),
                                  cl::NullRange, NULL, NULL);
  } catch (const cl::Error& e) {
    check_opencl_error("copy_triangular_transposed", e);
  }
}
}  // namespace math
}  // namespace stan

#endif
#endif
