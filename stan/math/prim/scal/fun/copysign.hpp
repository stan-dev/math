#ifndef STAN_MATH_PRIM_SCAL_FUN_COPYSIGN_HPP
#define STAN_MATH_PRIM_SCAL_FUN_COPYSIGN_HPP

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
 * @return magnitude of x with the sign of y
 */
template <class T, class U = T>
inline auto copysign(T const& t, U const& u) {
  return fabs(t) * sign(u);
}

}  // namespace math
}  // namespace stan
#endif
