#ifndef STAN_MATH_PRIM_MAT_FUN_SINH_HPP
#define STAN_MATH_PRIM_MAT_FUN_SINH_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap sinh() so that it can be vectorized.
 * @param x Angle in radians.
 * @tparam T Variable type.
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
 * @param x Container of variables.
 * @tparam T Container type.
 * @return Hyperbolic sine of each variable in x.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto sinh(const T& x) {
  return apply_scalar_unary<sinh_fun, T>::apply(x);
}

/**
 * Version of sinh() that accepts Eigen Matrix ar matrix expressions.
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Hyperbolic sine of each variable in x.
 */
template <typename Derived, typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto sinh(const Eigen::MatrixBase<Derived>& x){
  return x.derived().array().sinh().matrix();
}

/**
 * Version of acos() that accepts Eigen Array ar array expressions.
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Hyperbolic sine of each variable in x.
 */
template <typename Derived, typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto sinh(const Eigen::ArrayBase<Derived>& x){
  return x.derived().sinh();
}

}  // namespace math
}  // namespace stan

#endif
