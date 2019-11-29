#ifndef STAN_MATH_PRIM_MAT_FUN_ROUND_HPP
#define STAN_MATH_PRIM_MAT_FUN_ROUND_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/round.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap round() so it can be vectorized.
 * @param x Argument variable.
 * @tparam T Argument type.
 * @return Rounded value of x.
 */
struct round_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using stan::math::round;
    return round(x);
  }
};

/**
 * Vectorized version of round.
 * @param x Container.
 * @tparam T Container type.
 * @return Rounded value of each value in x.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto round(const T& x) {
  return apply_scalar_unary<round_fun, T>::apply(x);
}

/**
 * Version of round() that accepts Eigen Matrix ar matrix expressions.
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Rounded value of each value in x.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto round(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().round().matrix();
}

}  // namespace math
}  // namespace stan

#endif
