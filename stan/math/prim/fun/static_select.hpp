#ifndef STAN_MATH_PRIM_FUN_STATIC_SELECT_HPP
#define STAN_MATH_PRIM_FUN_STATIC_SELECT_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * If the condition is true returns the first value. Otherwise returns the
 * second value.
 * @tparam Condition the condition
 * @tparam T1 type of the first value
 * @tparam T2 type of the second value
 * @param a first value
 * @param b second value
 * @return first or second value, depending on the condition
 */
template <bool Condition, typename T1, typename T2,
          std::enable_if_t<Condition>* = nullptr>
T1 static_select(T1&& a, T2&& b) {
  return std::forward<T1>(a);
}

template <bool Condition, typename T2, typename T1,
          std::enable_if_t<!Condition>* = nullptr>
T2 static_select(T1&& a, T2&& b) {
  return std::forward<T2>(b);
}

}  // namespace math
}  // namespace stan

#endif
