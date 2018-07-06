#ifndef STAN_MATH_GPU_COPY_TRIANGULAR_OPENCL_HPP
#define STAN_MATH_GPU_COPY_TRIANGULAR_OPENCL_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/matrix_gpu.hpp>
#include <CL/cl.hpp>

namespace stan {
namespace math {
enum triangularity { LOWER = 0, UPPER = 1, NONE = 2 };

/**
 * Copies the lower or upper
 * triangular of the source matrix to
 * the destination matrix.
 * Both matrices are stored on the GPU.
 *
 * @param src the source matrix
 * @param lower_upper enum to describe
 * which part of the matrix to copy:
 * LOWER - copies the lower triangular
 * UPPER - copes the upper triangular
 *
 * @return the matrix with the copied content
 *
 */
inline matrix_gpu copy_triangular(matrix_gpu& src, triangularity lower_upper) {
  if (src.size() == 0) {
    return src;
  }
  if (src.size() == 1) {
    return src;
  }
  matrix_gpu dst(src.rows(), src.cols());
  cl::Kernel kernel = opencl_context.get_kernel("copy_triangular");
  cl::CommandQueue cmdQueue = opencl_context.queue();
  try {
    kernel.setArg(0, dst.buffer());
    kernel.setArg(1, src.buffer());
    kernel.setArg(2, dst.rows());
    kernel.setArg(3, dst.cols());
    kernel.setArg(4, lower_upper);
    cmdQueue.enqueueNDRangeKernel(kernel, cl::NullRange,
                                  cl::NDRange(dst.rows(), dst.cols()),
                                  cl::NullRange, NULL, NULL);
  } catch (const cl::Error& e) {
    check_opencl_error("copy_triangular", e);
  }
  return dst;
}
}  // namespace math
}  // namespace stan

#endif
#endif
