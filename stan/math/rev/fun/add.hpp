#ifndef STAN_MATH_REV_FUN_ADD_HPP
#define STAN_MATH_REV_FUN_ADD_HPP

#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Return the sum of the specified scalars.
 *
 * @tparam ScalarA type of the first scalar
 * @tparam ScalarB type of the second scalar
 * @param a first scalar
 * @param b second scalar
 * @return the sum of the scalars
 */
template <typename T1, typename T2, require_any_st_var<T1, T2>* = nullptr>
inline auto add(const T1& a, const T2& b) {
  return a + b;
}

}  // namespace math
}  // namespace stan

#endif
