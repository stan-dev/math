#ifndef STAN_MATH_PRIM_CORE_OPERATOR_DIVISION_HPP
#define STAN_MATH_PRIM_CORE_OPERATOR_DIVISION_HPP

#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/prim/meta/is_stan_scalar.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

namespace internal {
/** Evaluate a*b/c^2 without squaring c or losing a small factor before it
 * is multiplied by a large one. Branch on values before building AD nodes.
 */
template <typename T_a, typename T_b, typename T_c>
inline return_type_t<T_a, T_b, T_c> multiply_inv_square(const T_a& a,
                                                        const T_b& b,
                                                        const T_c& c) {
  const double av = std::fabs(value_of_rec(a));
  const double bv = std::fabs(value_of_rec(b));
  const double cv = std::fabs(value_of_rec(c));
  if (cv == 0 || !std::isfinite(cv)) {
    return (a * b) / (c * c);
  }
  if (av == 0 || bv == 0) {
    return ((a * b) / c) / c;
  }
  // Divide the larger factor first when shrinking, the smaller when growing.
  const bool a_first = (cv >= 1) == (av >= bv);
  const double product = a_first ? (av / cv) * bv : (bv / cv) * av;
  if (!std::isfinite(product) || (cv < 1 && !std::isnormal(product))) {
    return (a / c) * (b / c);
  }
  if (a_first) {
    return ((a / c) * b) / c;
  }
  return ((b / c) * a) / c;
}

/**
 * Return the quotient of the specified arguments.  At least one of the
 * arguments must be a complex number.
 *
 * @tparam U type of first argumentx
 * @tparam V type of second argument
 * @param[in] lhs first argument
 * @param[in] rhs second argument
 * @return quotient of the arguments
 */
template <typename U, typename V>
inline complex_return_t<U, V> complex_divide(const U& lhs, const V& rhs) {
  complex_return_t<U, V> y(lhs);
  y /= rhs;
  return y;
}
}  // namespace internal

/**
 * Return the quotient of the arguments.
 *
 * @tparam U value type of first argument
 * @tparam V value type of second argument
 * @param x first argument
 * @param y second argument
 * @return quotient of the arguments
 */
template <typename U, typename V, require_all_stan_scalar_t<U, V>* = nullptr>
inline complex_return_t<U, V> operator/(const std::complex<U>& x,
                                        const std::complex<V>& y) {
  return internal::complex_divide(x, y);
}

/**
 * Return the quotient of the arguments.
 *
 * @tparam U value type of first argument
 * @tparam V type of second argument
 * @param x first argument
 * @param y second argument
 * @return quotient of the arguments
 */
template <typename U, typename V, require_all_stan_scalar_t<U, V>* = nullptr>
inline complex_return_t<U, V> operator/(const std::complex<U>& x, const V& y) {
  return internal::complex_divide(x, y);
}

/**
 * Return the quotient of the arguments.
 *
 * @tparam U type of first argument
 * @tparam V value type of second argument
 * @param x first argument
 * @param y second argument
 * @return quotient of the arguments
 */
template <typename U, typename V, require_all_stan_scalar_t<U, V>* = nullptr>
inline complex_return_t<U, V> operator/(const U& x, const std::complex<V>& y) {
  return internal::complex_divide(x, y);
}

}  // namespace math
}  // namespace stan

#endif
