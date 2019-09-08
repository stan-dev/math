#ifndef STAN_MATH_PRIM_MAT_FUN_ADD_HPP
#define STAN_MATH_PRIM_MAT_FUN_ADD_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/err/check_matching_dims.hpp>
#include <boost/math/tools/promotion.hpp>

namespace stan {
namespace math {

/**
 * Return the sum of the specified matrices.  The two matrices
 * must have the same dimensions.
 * @tparam T1 Scalar type of first matrix.
 * @tparam T2 Scalar type of second matrix.
 * @tparam R Row type of matrices.
 * @tparam C Column type of matrices.
 * @param m1 First matrix.
 * @param m2 Second matrix.
 * @return Sum of the matrices.
 * @throw std::invalid_argument if m1 and m2 do not have the same
 * dimensions.
 */
template <typename Mat1, typename Mat2, require_all_eigen_t<is_stan_scalar, Mat1, Mat2>...>
inline auto add(Mat1&& m1, Mat2&& m2) {
  check_matching_dims("add", "m1", m1, "m2", m2);
  return (m1 + m2).eval();
}

/**
 * Return the sum of the specified matrix and specified scalar.
 *
 * @tparam T1 Type of matrix.
 * @tparam T2 Type of scalar.
 * @param m Matrix.
 * @param c Scalar.
 * @return The matrix plus the scalar.
 */
 template <typename T1, typename T2,
 require_eigen_t<is_stan_scalar, T1>...,
 require_not_eigen_t<is_stan_scalar, T2>...>
inline auto add(T1&& m, T2&& c) {
  return (m.array() + c).eval();
}

/**
 * Return the sum of the specified scalar and specified matrix.
 *
 * @tparam T1 Type of scalar.
 * @tparam T2 Type of matrix.
 * @param c Scalar.
 * @param m Matrix.
 * @return The scalar plus the matrix.
 */
 template <typename T1, typename T2,
 require_not_eigen_t<is_stan_scalar, T1>...,
 require_eigen_t<is_stan_scalar, T2>...>
inline auto add(T1&& c, T2& m) {
  return (c + m.array()).eval();
}

}  // namespace math
}  // namespace stan
#endif
