#ifndef STAN_MATH_PRIM_FUN_SQRT_HPP
#define STAN_MATH_PRIM_FUN_SQRT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the square root of the specified argument.  This
 * version is required to disambiguate <code>sqrt(int)</code>.
 *
 * @param[in] x Argument.
 * @return Square root of x.
 */
inline double sqrt(int x) { return std::sqrt(x); }

/**
 * Structure to wrap sqrt() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Square root of x.
 */
struct sqrt_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::sqrt;
    return sqrt(x);
  }
};

/**
 * Vectorized version of sqrt().
 *
 * @tparam T type of container
 * @param x container
 * @return Square root of each value in x.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto sqrt(const T& x) {
  return apply_scalar_unary<sqrt_fun, T>::apply(x);
}

/**
 * Version of sqrt() that accepts Eigen Matrix or matrix expressions.
 *
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Square root of each value in x.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto sqrt(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().sqrt().matrix().eval();
}

namespace internal {
/**
 * Return the square root of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return square root of the argument
 */
template <typename V>
inline std::complex<V> complex_sqrt(const std::complex<V>& z) {
  auto m = sqrt(hypot(z.real(), z.imag()));
  auto at = 0.5 * atan2(z.imag(), z.real());
  return {m * cos(at), m * sin(at)};
}
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
