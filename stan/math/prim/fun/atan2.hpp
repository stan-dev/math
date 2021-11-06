#ifndef STAN_MATH_PRIM_FUN_ATAN2_HPP
#define STAN_MATH_PRIM_FUN_ATAN2_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the principal value of the arg tangent of y/x,
 * expressed in radians.
 *
 * @tparam T1 type of first argument (must be arithmetic)
 * @tparam T2 type of second argument (must be arithmetic)
 * @param y first argument
 * @param x second argument
 * @return `atan2()` value of the arguments
 */
template <typename T1, typename T2, require_all_arithmetic_t<T1, T2>* = nullptr>
return_type_t<T1, T2> atan2(T1 y, T2 x) {
  return std::atan2(y, x);
}

}  // namespace math
}  // namespace stan

#endif
