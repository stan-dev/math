#ifndef STAN_MATH_PRIM_SCAL_FUN_COPYSIGN_HPP
#define STAN_MATH_PRIM_SCAL_FUN_COPYSIGN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the negation of the first argument if the first and second
 * argument have different signs, otherwise return a copy of the first
 * argument.  For the sake of this function, zero is considered
 * positive.  This function uses negation rather than literally copying
 * signs to preserve derivatives.
 *
 * Overload of `std::copysign` from `cmath` for argument-dependent
 * lookup.
 *
 * @tparam T type of first argument
 * @tparam U type of second argument
 * @param[in] x first complex argument
 * @param[in] y second complex argument
 * @return copy of the first argument, negated if necessary to match
 * the sign of the second argument
 */
template <typename T, typename U>
inline T copysign(const T& x, const U& y) {
  return (std::signbit(value_of_rec(x)) != std::signbit(value_of_rec(y))) ? -x
                                                                          : x;
}

/**
 * Return the negation of the first argument if the first and second
 * arguments have different signs and the first argument is not zero,
 * otherwise return a copy of the first argument.
 *
 * @tparam T type of first argument
 * @tparam U type of second argument
 * @param[in] x first complex argument
 * @param[in] y second complex argument
 * @return copy of the first argument, negated if necessary to match
 * the sign of the second argument
 * @see copysign
 */
template <typename T, typename U>
inline T copysign_non_zero(const T& x, const U& y) {
  return x != 0 ? stan::math::copysign(x, y) : x;
}

/**
 * Return the complex number composed of the real and complex parts
 * with signs copied from the real and complex parts of the first
 * arguments to the real and complex parts of the second.
 *
 * This is an overload of the standard libary `copysign` for complex
 * numbers that will be used with argument-dependent lookup.  Rather
 * than using the standard library `copysign`, it uses
 * `copysign_non_zero`, which does not change sign if the reference
 * value is zero (`-0.0` or `0.0`).
 *
 * @tparam T value type of first argument
 * @tparam U value type of second argument
 * @param[in] x first complex argument
 * @param[in] y second complex argument
 * @return copy of second argument, with components negated if
 * necessary to match sign of first argument
 */
template <typename T, typename U>
inline std::complex<T> copysign(const std::complex<T>& x,
                                const std::complex<U>& y) {
  return {copysign_non_zero(x.real(), y.real()),
          copysign_non_zero(x.imag(), y.imag())};
}

}  // namespace math
}  // namespace stan

#endif
