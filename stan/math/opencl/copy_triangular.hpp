#ifndef STAN_MATH_OPENCL_COPY_TRIANGULAR_HPP
#define STAN_MATH_OPENCL_COPY_TRIANGULAR_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernels/copy_triangular.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/prim/meta.hpp>
#include <cl.hpp>

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
 * matrix_cl_view::Lower - copies the lower triangular
 * matrix_cl_view::Upper - copes the upper triangular
 *
 * @return the matrix with the copied content
 *
 */
template <matrix_cl_view matrix_view = matrix_cl_view::Entire, typename T,
          typename = require_arithmetic_t<T>>
inline matrix_cl<T> copy_triangular(const matrix_cl<T>& src) {
  if (src.size() == 0 || src.size() == 1) {
    matrix_cl<T> dst(src);
    return dst;
  }
  matrix_cl_view dst_view = both(matrix_view, src.view());
  matrix_cl<T> dst(src.rows(), src.cols(), dst_view);
  try {
    opencl_kernels::copy_triangular(cl::NDRange(dst.rows(), dst.cols()), dst,
                                    src, dst.rows(), dst.cols(), dst_view);
  } catch (const cl::Error& e) {
    check_opencl_error("copy_triangular", e);
  }
  return dst;
}
}  // namespace math
}  // namespace stan

#endif
#endif
