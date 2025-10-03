#ifndef STAN_MATH_PRIM_META_IS_FLOATING_POINT_HPP
#define STAN_MATH_PRIM_META_IS_FLOATING_POINT_HPP

#include <type_traits>

namespace stan {

/**
 * Checks if decayed type is a floating point type
 * @tparam The type to check
 * @ingroup type_trait
 */
template <typename T>
using is_floating_point = std::is_floating_point<std::decay_t<T>>;

template <typename T>
inline constexpr bool is_floating_point_v = stan::is_floating_point<T>::value;

template <typename... Types>
inline constexpr bool is_all_floating_point_v
    = (stan::is_floating_point_v<Types> && ...);

template <typename... Types>
inline constexpr bool is_any_floating_point_v
    = (stan::is_floating_point_v<Types> || ...);

}  // namespace stan

#endif
