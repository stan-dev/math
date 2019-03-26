#ifndef STAN_MATH_PRIM_SCAL_FUN_SUM_HPP
#define STAN_MATH_PRIM_SCAL_FUN_SUM_HPP

namespace stan {
namespace math {

/**
 * Calculates sum of values given as input. For scalar inputs this is identity
 * function.
 *
 * @tparam T Type of element.
 * @param v Specified value.
 * @return Same value (the sum of one value).
 */
template <typename T>
inline T sum(const T& v) {
  return v;
}

}  // namespace math
}  // namespace stan
#endif
