#ifndef STAN_MATH_PRIM_SCAL_FUN_COPYSIGN_HPP
#define STAN_MATH_PRIM_SCAL_FUN_COPYSIGN_HPP

#include <stan/math/prim/scal/fun/signbit.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the magnitude of t with the sign of u
 *
 * Either T or U must be a stan type for ADL
 *
 * Needed for libc++'s implementation of
 * complex multiplication on stan types
 *
 * @tparam T type of object
 * @tparam U type of object
 * @param t magnitude reference
 * @param u sign reference
 * @return magnitude of t with the sign of u
 */
template <class T, class U = T>
inline auto copysign(T const& t, U const& u) {
  using std::fabs;
  using std::signbit;
  auto sign(signbit(u) ? -1 : 1);
  return fabs(t) * sign;
}

}  // namespace math
}  // namespace stan
#endif
