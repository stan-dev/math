#ifndef STAN_MATH_PRIM_FUN_SINH_HPP
#define STAN_MATH_PRIM_FUN_SINH_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap sinh() so that it can be vectorized.
 *
 * @tparam T type of argument
 * @param x angle in radians
 * @return Hyperbolic sine of x.
 */
struct sinh_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::sinh;
    return sinh(x);
  }
};

/**
 * Vectorized version of sinh().
 *
 * @tparam T type of container
 * @param x container
 * @return Hyperbolic sine of each variable in x.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto sinh(const T& x) {
  return apply_scalar_unary<sinh_fun, T>::apply(x);
}

/**
 * Version of sinh() that accepts Eigen Matrix or matrix expressions.
 *
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Hyperbolic sine of each variable in x.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto sinh(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().sinh().matrix().eval();
}

/**
 * Version of acos() that accepts Eigen Array or array expressions.
 *
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Hyperbolic sine of each variable in x.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto sinh(const Eigen::ArrayBase<Derived>& x) {
  return x.derived().sinh().eval();
}

}  // namespace math
}  // namespace stan

#endif
