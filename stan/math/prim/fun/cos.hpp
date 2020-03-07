#ifndef STAN_MATH_PRIM_FUN_COS_HPP
#define STAN_MATH_PRIM_FUN_COS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap cos() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x angle in radians
 * @return Cosine of x.
 */
struct cos_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::cos;
    return cos(x);
  }
};

/**
 * Vectorized version of cos().
 *
 * @tparam T type of container
 * @param x angles in radians
 * @return Cosine of each value in x.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto cos(const T& x) {
  return apply_scalar_unary<cos_fun, T>::apply(x);
}

/**
 * Version of cos() that accepts Eigen Matrix or matrix expressions.
 *
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Cosine of each value in x.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto cos(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().cos().matrix().eval();
}

}  // namespace math
}  // namespace stan

#endif
