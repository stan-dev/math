#ifndef STAN_MATH_PRIM_SCAL_FUN_COPYSIGN_HPP
#define STAN_MATH_PRIM_SCAL_FUN_COPYSIGN_HPP

namespace stan {
namespace math {

/**
 * Returns the magnitude of x with the sign of y
 * 
 * Needed for libc++'s implementation of
 * complex multiplication on stan types
 *
 * @tparam T type of stan object
 * @param x magnitude reference
 * @param y sign reference
 * @return magnitude of x with the sign of y
 */
template <class T>
inline auto
copysign(T const& x, T const& y) {
  return fabs(x) * sign(y);
}

}  // namespace math
}  // namespace stan
#endif
