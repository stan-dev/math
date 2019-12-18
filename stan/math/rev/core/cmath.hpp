#ifndef STAN_MATH_REV_CORE_CMATH_HPP
#define STAN_MATH_REV_CORE_CMATH_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/meta.hpp>


namespace stan {
namespace math {

/**
 * Return `true` if the specified argument is negative and `false`
 * otherwise.
 *
 * Overloads `std::signbit` from `<cmath>` for argument-dependent
 * lookup.
 *
 * @tparam T type of argument
 * @param[in] v argument
 * @return `true` if the argument is negative
 */
template <typename T, require_autodiff_t<T>...>
inline bool signbit(T&& v) {
  using std::signbit;
  return signbit(v.val());
}

/**
 * Return true if specified argument is infinite (positive or
 * negative).
 *
 * Overloads `std::isinf` from `<cmath>` for argument-dependent
 * lookup.
 *
 * @tparam T type of argument
 * @param[in] v argument
 * @return true if argument is infinite
 */
template <typename T, require_autodiff_t<T>...>
inline bool isinf(T&& v) {
  using std::isinf;
  return isinf(v.val());
}

/**
 * Return true if specified argument is finite (not infinite and not
 * not-a-number).
 *
 * Overloads `std::isfinite` from `<cmath>` for argument-dependent
 * lookup.
 *
 * @tparam T type of argument
 * @param[in] v argument
 * @return true if argument is finite
 */
template <typename T, require_autodiff_t<T>...>
inline bool isfinite(T&& v) {
  using std::isfinite;
  return isfinite(v.val());
}

/**
 * Return true if specified argument is not-a-number.
 *
 * Overloads `std::isnan` from `<cmath>` for argument-dependent
 * lookup.
 *
 * @tparam T type of argument
 * @param[in] v argument
 * @return true if argument is not-a-number
 */
template <typename T, require_autodiff_t<T>...>
inline bool isnan(T&& v) {
  using std::isnan;
  return isnan(v.val());
}

/**
 * Return true if specified argument is normal.  A number is normal if
 * it is finite, non-zero and not subnormal.
 *
 * Overloads `std::isnormal` from `<cmath>` for argument-dependent
 * lookup.
 *
 * @tparam T type of argument
 * @param[in] v argument
 * @return true if argument is normal
 */
template <typename T, require_autodiff_t<T>...>
inline bool isnormal(T&& v) {
  using std::isnormal;
  return isnormal(v.val());
}

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
 * @return copy of second argument, negated if necessary to match sign
 * of first argument
 */
template <typename T, typename U>
inline T copysign(const T& x, const U& y) {
  // 0 is considered positive
  return (x < 0 && y >= 0) || (x > 0 && y < 0) ? -x : x;
}
}  // namespace math
}  // namespace stan

#endif
