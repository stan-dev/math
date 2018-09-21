#ifndef STAN_MATH_GPU_TRANSPOSE_HPP
#define STAN_MATH_GPU_TRANSPOSE_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/matrix_gpu.hpp>
#include <stan/math/gpu/kernels/transpose.hpp>
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
  cl::CommandQueue cmdQueue = opencl_context.queue();
  try {
    opencl_kernels::transpose(cl::NDRange(src.rows(), src.cols()), dst.buffer(),
                              src.buffer(), src.rows(), src.cols());
  } catch (const cl::Error& e) {
    check_opencl_error("transpose", e);
  }
  return dst;
}
}  // namespace math
}  // namespace stan

#endif
#endif
