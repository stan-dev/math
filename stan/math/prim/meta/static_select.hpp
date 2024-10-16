#ifndef STAN_MATH_PRIM_META_STATIC_SELECT_HPP
#define STAN_MATH_PRIM_META_STATIC_SELECT_HPP

#include <type_traits>
#include <utility>

namespace stan {
namespace math {

/**
 * Returns one of the arguments that can be of different type, depending on the
 * compile time condition. Can often be used in place of c++17 `if constexpr`.
 *
 * @tparam Condition condition
 * @tparam T1 type of the first argument
 * @tparam T2 type of the second argument
 *
 * @param a first argument; returned if the `Condition` is true
 * @param b second argument; returned if the `Condition` is false
 * @return One of the arguemtns, depending on the `Condition`
 */
template <bool Condition, typename T1, typename T2,
          require_t<std::bool_constant<Condition>>* = nullptr>
inline T1 static_select(T1&& a, T2&& b) {
  return std::forward<T1>(a);
}

template <bool Condition, typename T1, typename T2,
          require_not_t<std::bool_constant<Condition>>* = nullptr>
inline T2 static_select(T1&& a, T2&& b) {
  return std::forward<T2>(b);
}

}  // namespace math
}  // namespace stan

#endif  // STATIC_SELECT_HPP
