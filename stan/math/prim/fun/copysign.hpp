#ifndef STAN_MATH_PRIM_SCAL_FUN_COPYSIGN_HPP
#define STAN_MATH_PRIM_SCAL_FUN_COPYSIGN_HPP

#include <stan/math/prim/meta.hpp>

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
 * @tparam ADType type of first argument
 * @tparam U type of second argument
 * @param[in] x first complex argument
 * @param[in] y second complex argument
 * @return copy of second argument, negated if necessary to match sign
 * of first argument
 */
template <typename ADType, typename U>
inline ADType copysign(const ADType& x, const U& y) {
  // +0 is positive; second condition handles -0
  return (x < 0 && y >= 0) || (x >= 0 && y < 0) ? -x : x;
}
}  // namespace math
}  // namespace stan

#endif
