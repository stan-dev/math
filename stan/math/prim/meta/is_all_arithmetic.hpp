#ifndef STAN_MATH_PRIM_META_IS_ALL_ARITHMETIC_HPP
#define STAN_MATH_PRIM_META_IS_ALL_ARITHMETIC_HPP

#include <stan/math/prim/meta/scalar_type.hpp>
#include <type_traits>
#include <tuple>

namespace stan {
template <typename T>
using is_arithmetic = std::is_arithmetic<scalar_type_t<T>>;

template <typename T>
inline constexpr bool is_arithmetic_v = is_arithmetic<std::decay_t<T>>::value;

namespace internal {

template <typename... Types>
struct is_all_arithmetic_scalar_impl
    : std::conjunction<is_arithmetic<std::decay_t<Types>>...> {};

template <typename... Types>
struct is_all_arithmetic_scalar_impl<std::tuple<Types...>>
    : std::conjunction<is_all_arithmetic_scalar_impl<std::decay_t<Types>>...> {
};
template <typename T, typename... VecArgs>
struct is_all_arithmetic_scalar_impl<std::vector<T, VecArgs...>>
    : std::conjunction<is_all_arithmetic_scalar_impl<std::decay_t<T>>> {};
}  // namespace internal

template <typename... Types>
struct is_all_arithmetic_scalar
    : std::conjunction<
          internal::is_all_arithmetic_scalar_impl<std::decay_t<Types>>...> {};

template <typename... Types>
inline constexpr bool is_all_arithmetic_scalar_v
    = is_all_arithmetic_scalar<std::decay_t<Types>...>::value;

}  // namespace stan

#endif
