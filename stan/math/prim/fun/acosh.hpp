#ifndef STAN_MATH_PRIM_FUN_ACOSH_HPP
#define STAN_MATH_PRIM_FUN_ACOSH_HPP

#include <stan/math/prim/fun/is_nan.hpp>
#include <stan/math/prim/fun/is_inf.hpp>
#include <stan/math/prim/err/check_greater_or_equal.hpp>
#include <cmath>
#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>










namespace stan {
namespace math {

/**
 * Return the inverse hyperbolic cosine of the specified value.
 * Returns nan for nan argument.
 *
 * @param[in] x Argument.
 * @return Inverse hyperbolic cosine of the argument.
 * @throw std::domain_error If argument is less than 1.
 */
inline double acosh(double x) {
  if (is_nan(x)) {
    return x;
  } else {
    check_greater_or_equal("acosh", "x", x, 1.0);
#ifdef _WIN32
    if (is_inf(x))
      return x;
#endif
    return std::acosh(x);
  }
}

/**
 * Integer version of acosh.
 *
 * @param[in] x Argument.
 * @return Inverse hyperbolic cosine of the argument.
 * @throw std::domain_error If argument is less than 1.
 */
inline double acosh(int x) {
  if (is_nan(x)) {
    return x;
  } else {
    check_greater_or_equal("acosh", "x", x, 1);
#ifdef _WIN32
    if (is_inf(x))
      return x;
#endif
    return std::acosh(x);
  }
}













/**
 * Structure to wrap acosh() so it can be vectorized.
 */
struct acosh_fun {
  /**
   * Return the inverse hypberbolic cosine of the specified argument.
   *
   * @param x Argument.
   * @return Inverse hyperbolic cosine of the argument.
   * @tparam T Argument type.
   */
  template <typename T>
  static inline T fun(const T& x) {
    return acosh(x);
  }
};

/**
 * Return the elementwise application of <code>acosh()</code> to
 * specified argument container.  The return type promotes the
 * underlying scalar argument type to double if it is an integer,
 * and otherwise is the argument type.
 *
 * @tparam T Container type.
 * @param x Container.
 * @return Elementwise acosh of members of container.
 */
template <typename T>
inline typename apply_scalar_unary<acosh_fun, T>::return_t acosh(const T& x) {
  return apply_scalar_unary<acosh_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
