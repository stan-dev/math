#ifndef STAN_MATH_GPU_COPY_TRIANGULAR_HPP
#define STAN_MATH_GPU_COPY_TRIANGULAR_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/constants.hpp>
#include <stan/math/gpu/matrix_gpu.hpp>
#include <CL/cl.hpp>

namespace stan {
namespace math {

/**
 * Copies the lower or upper
 * triangular of the source matrix to
 * the destination matrix.
 * Both matrices are stored on the GPU.
 *
 * @param src the source matrix
 * @tparam triangular_map int to describe
 * which part of the matrix to copy:
 * Lower - copies the lower triangular
 * Upper - copes the upper triangular
 *
 * @return the matrix with the copied content
 *
 */
template <int triangular_map = gpu::Lower>
inline matrix_gpu copy_triangular(const matrix_gpu& src) {
  if (src.size() == 0 || src.size() == 1) {
    matrix_gpu dst(src);
    return dst;
  }
  matrix_gpu dst(src.rows(), src.cols());
  cl::CommandQueue cmdQueue = opencl_context.queue();
  try {
    kernel_cl kernel("copy_triangular");
    kernel.set_args(dst.buffer(), src.buffer(), dst.rows(), dst.cols(),
                    triangular_map);
    cmdQueue.enqueueNDRangeKernel(kernel.compiled_, cl::NullRange,
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
