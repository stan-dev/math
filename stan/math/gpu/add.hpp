#ifndef STAN_MATH_GPU_ADD_HPP
#define STAN_MATH_GPU_ADD_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/matrix_gpu.hpp>
#include <stan/math/gpu/kernels/add.hpp>
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
  cl::CommandQueue cmdQueue = opencl_context.queue();
  try {
    opencl_kernels::add(cl::NDRange(A.rows(), A.cols()), C.buffer(), A.buffer(),
                        B.buffer(), A.rows(), A.cols());
  } catch (const cl::Error& e) {
    check_opencl_error("add", e);
  }
  return C;
}

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
inline auto operator+(const matrix_gpu& A, const matrix_gpu& B) {
  return add(A, B);
}
}  // namespace math
}  // namespace stan

#endif
#endif
