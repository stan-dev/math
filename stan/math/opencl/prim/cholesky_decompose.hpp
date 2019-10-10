#ifndef STAN_MATH_OPENCL_PRIM_CHOLESKY_DECOMPOSE_HPP
#define STAN_MATH_OPENCL_PRIM_CHOLESKY_DECOMPOSE_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/cholesky_decompose.hpp>
#include <stan/math/opencl/copy_triangular.hpp>
#include <stan/math/prim/meta.hpp>
#include <cl.hpp>
#include <algorithm>
#include <cmath>

namespace stan {
namespace math {
/**
 * Returns the lower-triangular Cholesky factor (i.e., matrix
 * square root) of the specified square, symmetric matrix on the OpenCL device.
 * The return value \f$L\f$ will be a lower-traingular matrix such that the
 * original matrix \f$A\f$ is given by <p>\f$A = L \times L^T\f$.
 * @param A Input square matrix
 * @return Square root of matrix.
 * @throw std::domain_error if m is not a symmetric matrix or
 *   if m is not positive definite (if m has more than 0 elements)
 */
template <typename T, typename = require_floating_point_t<T>>
inline matrix_cl<T> cholesky_decompose(matrix_cl<T>& A) {
  check_square("cholesky_decompose", "A", A);
  check_symmetric("cholesky_decompose", "A", A);
  matrix_cl<T> res = copy_cl(A);
  if (res.rows() == 0) {
    return res;
  }
  opencl::cholesky_decompose(res);
  // check_pos_definite on matrix_cl is check_nan + check_diagonal_zeros
  check_nan("cholesky_decompose (OpenCL)", "A", res);
  check_diagonal_zeros("cholesky_decompose (OpenCL)", "A", res);
  res.template zeros_strict_tri<matrix_cl_view::Upper>();
  return res;
}
}  // namespace math
}  // namespace stan

#endif
#endif
