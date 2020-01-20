#ifndef STAN_MATH_PRIM_FUN_SUBTRACT_HPP
#define STAN_MATH_PRIM_FUN_SUBTRACT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the result of subtracting the second specified matrix
 * from the first specified matrix.
 *
 * @tparam T1 type of the first matrix or expression
 * @tparam T2 type of the second matrix or expression
 *
 * @param m1 First matrix or expression.
 * @param m2 Second matrix or expression.
 * @return Difference between first matrix and second matrix.
 */
template <typename T1, typename T2, typename = require_all_eigen_t<T1, T2>>
inline auto subtract(const T1& m1, const T2& m2) {
  check_matching_dims("subtract", "m1", m1, "m2", m2);
  return m1 - m2;
}

/**
 * Return the result of subtracting the specified matrix from the specified
 * matrix.
 *
 * @tparam T1 type of the scalar
 * @tparam T2 type of the matrix or expression
 * @param c Scalar.
 * @param m Matrix or expression.
 * @return The matrix minus the scalar.
 */
template <typename T1, typename T2, typename = require_stan_scalar_t<T1>,
          typename = require_eigen_t<T2>>
inline auto subtract(const T1 c, const T2& m) {
  return (c - m.array()).matrix();
}

/**
 * Return the result of subtracting the specified scalar from the specified
 * matrix.
 *
 * @tparam T1 type of the matrix or expression
 * @tparam T2 type of the scalar
 * @param m Matrix or expression.
 * @param c Scalar.
 * @return The matrix minus the scalar.
 */
template <typename T1, typename T2, typename = require_eigen_t<T1>,
          typename = require_stan_scalar_t<T2>>
inline auto subtract(const T1& m, const T2 c) {
  return (m.array() - c).matrix();
}

}  // namespace math
}  // namespace stan

#endif
