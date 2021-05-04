#ifndef STAN_MATH_PRIM_CORE_OPERATOR_NOT_EQUAL_HPP
#define STAN_MATH_PRIM_CORE_OPERATOR_NOT_EQUAL_HPP

#include <stan/math/prim/meta.hpp>
#include <complex>

namespace stan {
namespace math {

/**
 * Return `true` if the complex numbers have unequal imaginary or
 * complex parts.
 *
 * @tparam U value type of first argument
 * @tparam V value type of second argument
 * @param x first argument
 * @param y second argument
 * @return `true` if the arguments are equal
 */
template <typename U, typename V, typename = require_any_autodiff_t<U, V>>
inline bool operator!=(const std::complex<U>& x, const std::complex<V>& y) {
  return !(x.real() == y.real() && x.imag() == y.imag());
}

/**
 * Return `true` if the first argument's real part is unequal to the
 * second argument or the first argument's imaginary part is unequal
 * to zero.
 *
 * @tparam U value type of first argument
 * @tparam V type of second argument
 * @param x first argument
 * @param y second argument
 * @return `true` if the arguments are equal
 */
template <typename U, typename V, typename = require_any_autodiff_t<U, V>>
inline bool operator!=(const std::complex<U>& x, const V& y) {
  return !(x.real() == y && x.imag() == 0);
}

/**
 * Return `true` if the first argument is unequal to the real part of
 * the second argument or the imaginary part of the second argument
 * is nonzero.
 *
 * @tparam U type of first argument
 * @tparam V value type of second argument
 * @param x first argument
 * @param y second argument
 * @return `true` if the arguments are not equal
 */
template <typename U, typename V, typename = require_any_autodiff_t<U, V>>
inline bool operator!=(const U& x, const std::complex<V>& y) {
  return !(x == y.real() && 0 == y.imag());
}

}  // namespace math
}  // namespace stan

#endif
