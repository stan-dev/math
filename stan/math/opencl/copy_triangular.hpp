#ifndef STAN_MATH_OPENCL_COPY_TRIANGULAR_HPP
#define STAN_MATH_OPENCL_COPY_TRIANGULAR_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/constants.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernels/copy_triangular.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/prim/meta.hpp>
#include <CL/cl.hpp>

namespace stan {
namespace math {

/**
 * Copies the lower or upper
 * triangular of the source matrix to
 * the destination matrix.
 * Both matrices are stored on the OpenCL device.
 *
 * @param src the source matrix
 * @tparam triangular_map int to describe
 * which part of the matrix to copy:
 * TriangularViewCL::Lower - copies the lower triangular
 * TriangularViewCL::Upper - copes the upper triangular
 *
 * @return the matrix with the copied content
 *
 */
template <TriangularViewCL triangular_view = TriangularViewCL::Entire,
          typename T, typename = enable_if_arithmetic<T>>
inline matrix_cl<T> copy_triangular(const matrix_cl<T>& src) {
  if (src.size() == 0 || src.size() == 1) {
    matrix_cl<T> dst(src);
    return dst;
  }
  matrix_cl<T> dst(src.rows(), src.cols());
  try {
    opencl_kernels::copy_triangular(cl::NDRange(dst.rows(), dst.cols()), dst,
                                    src, dst.rows(), dst.cols(),
                                    triangular_view);
  } catch (const cl::Error& e) {
    check_opencl_error("copy_triangular", e);
  }
  return dst;
}
}  // namespace math
}  // namespace stan

#endif
#endif
