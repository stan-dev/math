#ifndef STAN_MATH_PRIM_META_IS_INTEGRAL_HPP
#define STAN_MATH_PRIM_META_IS_INTEGRAL_HPP

#include <type_traits>

namespace stan {

/**
 * Checks if decayed type is integral
 * @tparam The type to check
 * @ingroup type_trait
 */
template <typename T>
using is_integral = std::is_integral<std::decay_t<T>>;

template <typename T>
inline constexpr bool is_integral_v = stan::is_integral<T>::value;

template <typename... Types>
inline constexpr bool is_all_integral_v = (stan::is_integral_v<Types> && ...);

template <typename... Types>
inline constexpr bool is_any_integral_v = (stan::is_integral_v<Types> || ...);

}  // namespace stan

#endif
