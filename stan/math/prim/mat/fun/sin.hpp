#ifndef STAN_MATH_PRIM_MAT_FUN_SIN_HPP
#define STAN_MATH_PRIM_MAT_FUN_SIN_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap sin() so it can be vectorized.
 * @param x Angle in radians.
 * @tparam T Argument type.
 * @return Sine of x.
 */
struct sin_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::sin;
    return sin(x);
  }
};

/**
 * Vectorized version of sin().
 * @param x Container of angles in radians.
 * @tparam T Container type.
 * @return Sine of each value in x.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto sin(const T& x) {
  return apply_scalar_unary<sin_fun, T>::apply(x);
}

/**
 * Version of sin() that accepts Eigen Matrix ar matrix expressions.
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Sine of each value in x.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto sin(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().sin().matrix();
}

}  // namespace math
}  // namespace stan

#endif
