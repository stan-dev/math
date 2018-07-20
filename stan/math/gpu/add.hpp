#ifndef STAN_MATH_GPU_ADD_HPP
#define STAN_MATH_GPU_ADD_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/matrix_gpu.hpp>
#include <stan/math/gpu/err/check_matching_dims.hpp>
#include <CL/cl.hpp>

namespace stan {
namespace math {

/**
 * Matrix addition on the GPU
 *
 * @param A first matrix
 * @param B second matrix
 *
 * @return sum of A and B
 *
 * @throw <code>std::invalid_argument</code> if the
 * input matrices do not have matching dimensions
 *
 */
inline matrix_gpu add(const matrix_gpu& A, const matrix_gpu& B) {
  check_matching_dims("add", "A", A, "B", B);
  matrix_gpu C(A.rows(), A.cols());
  if (C.size() == 0) {
    return C;
  }
  cl::Kernel kernel = opencl_context.get_kernel("add");
  cl::CommandQueue cmdQueue = opencl_context.queue();
  try {
    opencl_context.set_kernel_args(kernel, C.buffer(), A.buffer(), B.buffer(),
                                   A.rows(), A.cols());
    cmdQueue.enqueueNDRangeKernel(kernel, cl::NullRange,
                                  cl::NDRange(A.rows(), A.cols()),
                                  cl::NullRange, NULL, NULL);
  } catch (const cl::Error& e) {
    check_opencl_error("add", e);
  }
  return C;
}
}  // namespace math
}  // namespace stan

#endif
#endif
