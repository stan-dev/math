#ifndef STAN_MATH_PRIM_FUN_LOG_HPP
#define STAN_MATH_PRIM_FUN_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the natural log of the specified argument.  This version
 * is required to disambiguate <code>log(int)</code>.
 *
 * @param[in] x Argument.
 * @return Natural log of argument.
 */
inline double log(int x) { return std::log(x); }

/**
 * Structure to wrap log() so that it can be vectorized.
 */
struct log_fun {
  /**
   * Return natural log of specified argument.
   *
   * @tparam T type of argument
   * @param[in] x argument
   * @return Natural log of x.
   */
  template <typename T>
  static inline T fun(const T& x) {
    using std::log;
    return log(x);
  }
};

/**
 * Return the elementwise natural log of the specified argument,
 * which may be a scalar or any Stan container of numeric scalars.
 * The return type is the same as the argument type.
 *
 * @tparam T type of container
 * @param[in] x container
 * @return Elementwise application of natural log to the argument.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto log(const T& x) {
  return apply_scalar_unary<log_fun, T>::apply(x);
}

/**
 * Version of log() that accepts Eigen Matrix or matrix expressions.
 *
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Arc cosine of each variable in the container, in radians.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto log(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().log().matrix().eval();
}

}  // namespace math
}  // namespace stan
#endif
