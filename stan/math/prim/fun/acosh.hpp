#ifndef STAN_MATH_PRIM_FUN_ACOSH_HPP
#define STAN_MATH_PRIM_FUN_ACOSH_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/copysign.hpp>
#include <stan/math/prim/fun/is_inf.hpp>
#include <stan/math/prim/fun/is_nan.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>
#include <complex>

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
   * Return the inverse hyperbolic cosine of the specified argument.
   *
   * @tparam T type of argument
   * @param x argument
   * @return Inverse hyperbolic cosine of the argument.
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
 * @tparam T type of container
 * @param x container
 * @return Elementwise acosh of members of container.
 */
template <
    typename T, require_not_var_matrix_t<T>* = nullptr,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto acosh(const T& x) {
  return apply_scalar_unary<acosh_fun, T>::apply(x);
}

namespace internal {

/**
 * Return the hyperbolic arc cosine of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return hyperbolic arc cosine of the argument
 */
template <typename V>
inline std::complex<V> complex_acosh(const std::complex<V>& z) {
  std::complex<double> y_d = acosh(value_of_rec(z));
  auto y = log(z + sqrt(z * z - 1));
  return copysign(y, y_d);
}

}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
