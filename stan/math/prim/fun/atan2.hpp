#ifndef STAN_MATH_PRIM_FUN_ATAN2_HPP
#define STAN_MATH_PRIM_FUN_ATAN2_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Computes the arc tangent of y/x using the signs of
 * arguments to determine the correct quadrant.
 *
 * @tparam T1 type of first argument (must be arithmetic)
 * @tparam T2 type of second argument (must be arithmetic)
 * @param y first argument
 * @param x second argument
 * @return arctangent of `y / x`
 */
template <typename T1, typename T2, require_all_arithmetic_t<T1, T2>* = nullptr>
double atan2(T1 y, T2 x) {
  return std::atan2(y, x);
}

}  // namespace math
}  // namespace stan

#endif
