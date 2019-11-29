#ifndef STAN_MATH_PRIM_MAT_FUN_TAN_HPP
#define STAN_MATH_PRIM_MAT_FUN_TAN_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap tan() so that it can be vectorized.
 * @param x Angle in radians.
 * @tparam T Variable type.
 * @return Tangent of x.
 */
struct tan_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::tan;
    return tan(x);
  }
};

/**
 * Vectorized version of tan().
 * @param x Container of angles in radians.
 * @tparam T Container type.
 * @return Tangent of each value in x.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto tan(const T& x) {
  return apply_scalar_unary<tan_fun, T>::apply(x);
}

/**
 * Version of tan() that accepts Eigen Matrix ar matrix expressions.
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Tangent of each value in x.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto tan(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().tan().matrix();
}

}  // namespace math
}  // namespace stan

#endif
