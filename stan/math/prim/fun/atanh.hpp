#ifndef STAN_MATH_PRIM_FUN_ATANH_HPP
#define STAN_MATH_PRIM_FUN_ATANH_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/copysign.hpp>
#include <stan/math/prim/fun/is_nan.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the inverse hyperbolic tangent of the specified value.
 * An argument of -1 returns negative infinity and an argument of 1
 * returns infinity.
 * Returns nan for nan argument.
 *
 * @param[in] x `Arithmetic` Argument.
 * @return Inverse hyperbolic tangent of the argument.
 * @throw std::domain_error If argument is not in [-1, 1].
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline double atanh(const T x) {
  if (is_nan(x)) {
    return x;
  } else {
    check_bounded("atanh", "x", x, -1.0, 1.0);
    return std::atanh(x);
  }
}

/**
 * Return the inverse hyperbolic tangent of the specified value.
 * An argument of -1 returns negative infinity and an argument of 1
 * returns infinity.
 * Returns nan for nan argument.
 *
 * @param[in] x `complex<Arithmetic>` Argument.
 * @return Inverse hyperbolic tangent of the argument.
 * @throw std::domain_error If argument is not in [-1, 1].
 */
template <typename T, require_complex_bt<std::is_arithmetic, T>* = nullptr>
inline auto atanh(const T x) {
  return std::atanh(x);
}

/**
 * Structure to wrap atanh() so it can be vectorized.
 */
struct atanh_fun {
  /**
   * Return the inverse hyperbolic tangent of the specified argument.
   *
   * @tparam T type of argument
   * @param x argument
   * @return Inverse hyperbolic tangent of the argument.
   */
  template <typename T>
  static inline auto fun(const T& x) {
    return atanh(x);
  }
};

/**
 * Return the elementwise application of <code>atanh()</code> to
 * specified argument container.  The return type promotes the
 * underlying scalar argument type to double if it is an integer,
 * and otherwise is the argument type.
 *
 * @tparam T type of container
 * @param x container
 * @return Elementwise atanh of members of container.
 */
template <typename T, require_ad_container_t<T>* = nullptr>
inline auto atanh(const T& x) {
  return apply_scalar_unary<atanh_fun, T>::apply(x);
}

/**
 * Return the elementwise application of <code>atanh()</code> to
 * specified argument container.  The return type promotes the
 * underlying scalar argument type to double if it is an integer,
 * and otherwise is the argument type.
 *
 * @tparam T type of container
 * @param x container
 * @return Elementwise atanh of members of container.
 */
template <typename Container,
          require_container_bt<std::is_arithmetic, Container>* = nullptr>
inline auto atanh(const Container& x) {
  return apply_scalar_unary<atanh_fun, Container>::apply(x);
}

namespace internal {
/**
 * Return the hyperbolic arc tangent of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return hyperbolic arc tangent of the argument
 */
template <typename V>
inline std::complex<V> complex_atanh(const std::complex<V>& z) {
  std::complex<double> y_d = atanh(value_of_rec(z));
  V one(1);
  auto y = 0.5 * (log(one + z) - log(one - z));
  return copysign(y, y_d);
}
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
