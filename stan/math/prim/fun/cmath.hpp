#ifndef STAN_MATH_PRIM_SCAL_FUN_CMATH_HPP
#define STAN_MATH_PRIM_SCAL_FUN_CMATH_HPP

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
 * @tparam ADType type of argument
 * @param[in] v argument
 * @return `true` if the argument is negative
 */
template <typename ADType, require_autodiff_t<ADType>...>
inline bool signbit(ADType&& v) {
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
 * @tparam ADType type of argument
 * @param[in] v argument
 * @return true if argument is infinite
 */
template <typename ADType, require_autodiff_t<ADType>...>
inline bool isinf(ADType&& v) {
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
 * @tparam ADType type of argument
 * @param[in] v argument
 * @return true if argument is finite
 */
template <typename ADType, require_autodiff_t<ADType>...>
inline bool isfinite(ADType&& v) {
  using std::isfinite;
  return isfinite(v.val());
}

/**
 * Return true if specified argument is not-a-number.
 *
 * Overloads `std::isnan` from `<cmath>` for argument-dependent
 * lookup.
 *
 * @tparam ADType type of argument
 * @param[in] v argument
 * @return true if argument is not-a-number
 */
template <typename ADType, require_autodiff_t<ADType>...>
inline bool isnan(ADType&& v) {
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
 * @tparam ADType type of argument
 * @param[in] v argument
 * @return true if argument is normal
 */
template <typename ADType, require_autodiff_t<ADType>...>
inline bool isnormal(ADType&& v) {
  using std::isnormal;
  return isnormal(v.val());
}

/**
 * Return the negation of the first argument if the first and second
 * argument have different signs, otherwise return a copy of the first
 * argument.  For the sake of this function, zero is considered
 * positive.  ADTypehis function uses negation rather than literally copying
 * signs to preserve derivatives.
 *
 * Overload of `std::copysign` from `cmath` for argument-dependent
 * lookup.
 *
 * @tparam ADType type of first argument
 * @tparam U type of second argument
 * @param[in] x first complex argument
 * @param[in] y second complex argument
 * @return copy of second argument, negated if necessary to match sign
 * of first argument
 */
template <typename ADType, typename U>
inline ADType copysign(const ADType& x, const U& y) {
  // 0 is considered positive
  return (x < 0 && y >= 0) || (x > 0 && y < 0) ? -x : x;
}
}  // namespace math
}  // namespace stan

#endif
