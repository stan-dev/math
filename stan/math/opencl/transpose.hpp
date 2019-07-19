#ifndef STAN_MATH_OPENCL_TRANSPOSE_HPP
#define STAN_MATH_OPENCL_TRANSPOSE_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernels/transpose.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>

#include <CL/cl.hpp>

namespace stan {
namespace math {
/**
 * Takes the transpose of the matrix on the OpenCL device.
 *
 * @param src the input matrix
 *
 * @return transposed input matrix
 *
 */
template <typename T,
          typename std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline matrix_cl<T> transpose(const matrix_cl<T>& src) {
  matrix_cl<T> dst(src.cols(), src.rows());
  if (dst.size() == 0)
    return dst;
  try {
    opencl_kernels::transpose(cl::NDRange(src.rows(), src.cols()), dst, src,
                              src.rows(), src.cols());
  } catch (const cl::Error& e) {
    check_opencl_error("transpose", e);
  }
  return dst;
}
}  // namespace math
}  // namespace stan

#endif
#endif
