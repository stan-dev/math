#ifndef STAN_MATH_PRIM_FUN_ROUND_HPP
#define STAN_MATH_PRIM_FUN_ROUND_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the closest integer to the specified argument, with
 * halfway cases rounded away from zero.
 *
 * @param x Argument.
 * @return The rounded value of the argument.
 */
inline double round(double x) { return std::round(x); }

/**
 * Return the closest integer to the specified argument, with
 * halfway cases rounded away from zero.
 *
 * @param x Argument.
 * @return The rounded value of the argument.
 */
inline double round(int x) { return std::round(x); }

/**
 * Structure to wrap round() so it can be vectorized.
 *
 * @tparam T type of argument
 * @param x argument variable
 * @return Rounded value of x.
 */
struct round_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return round(x);
  }
};

/**
 * Vectorized version of round.
 *
 * @tparam T type of container
 * @param x container
 * @return Rounded value of each value in x.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto round(const T& x) {
  return apply_scalar_unary<round_fun, T>::apply(x);
}

/**
 * Version of round() that accepts Eigen Matrix or matrix expressions.
 *
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Rounded value of each value in x.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto round(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().round().matrix().eval();
}

}  // namespace math
}  // namespace stan

#endif
