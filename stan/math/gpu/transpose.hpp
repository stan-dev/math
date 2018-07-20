#ifndef STAN_MATH_GPU_TRANSPOSE_HPP
#define STAN_MATH_GPU_TRANSPOSE_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/matrix_gpu.hpp>
#include <CL/cl.hpp>

namespace stan {
namespace math {
/**
 * Takes the transpose of the matrix on the GPU.
 *
 * @param src the input matrix
 *
 * @return transposed input matrix
 *
 */
inline matrix_gpu transpose(const matrix_gpu& src) {
  matrix_gpu dst(src.cols(), src.rows());
  if (dst.size() == 0)
    return dst;
  cl::Kernel kernel = opencl_context.get_kernel("transpose");
  cl::CommandQueue cmdQueue = opencl_context.queue();
  try {
    opencl_context.set_kernel_args(kernel, dst.buffer(), src.buffer(),
                                   src.rows(), src.cols());
    cmdQueue.enqueueNDRangeKernel(kernel, cl::NullRange,
                                  cl::NDRange(src.rows(), src.cols()),
                                  cl::NullRange, NULL, NULL);
  } catch (const cl::Error& e) {
    check_opencl_error("transpose", e);
  }
  return dst;
}
}  // namespace math
}  // namespace stan

#endif
#endif
