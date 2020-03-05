#ifndef STAN_MATH_PRIM_FUN_LOG10_HPP
#define STAN_MATH_PRIM_FUN_LOG10_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Structure to wrap log10() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Log base-10 of x.
 */
struct log10_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::log10;
    return log10(x);
  }
};

/**
 * Vectorized version of log10().
 *
 * @tparam T type of container
 * @param x container
 * @return Log base-10 applied to each value in x.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto log10(const T& x) {
  return apply_scalar_unary<log10_fun, T>::apply(x);
}

/**
 * Version of log10() that accepts Eigen Matrix or matrix expressions.
 *
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Log base-10 applied to each value in x.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto log10(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().log10().matrix().eval();
}

namespace internal {
/**
 * Return the base 10 logarithm of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return base 10 logarithm of the argument
 */
template <typename V>
inline std::complex<V> complex_log10(const std::complex<V>& z) {
  static const double inv_log_10 = 1 / std::log(10);
  return log(z) * inv_log_10;
}
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
