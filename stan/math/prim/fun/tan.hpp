#ifndef STAN_MATH_PRIM_FUN_TAN_HPP
#define STAN_MATH_PRIM_FUN_TAN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap tan() so that it can be vectorized.
 *
 * @tparam T type of argument
 * @param x angle in radians
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
 *
 * @tparam T type of container
 * @param x angles in radians
 * @return Tangent of each value in x.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto tan(const T& x) {
  return apply_scalar_unary<tan_fun, T>::apply(x);
}

/**
 * Version of tan() that accepts Eigen Matrix or matrix expressions.
 *
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Tangent of each value in x.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto tan(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().tan().matrix().eval();
}

}  // namespace math
}  // namespace stan

#endif
