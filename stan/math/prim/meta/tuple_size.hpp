#ifndef STAN_MATH_PRIM_META_TUPLE_SIZE_HPP
#define STAN_MATH_PRIM_META_TUPLE_SIZE_HPP

#include <stan/math/prim/meta/is_tuple.hpp>
#include <type_traits>
#include <cstddef>
#include <tuple>

namespace stan {

/**
 * Equivalent to std::tuple_size but returns 0 T is not a tuple
 * @tparam T type to get tuple size of
 * @ingroup type_trait
 */
template <typename T, typename = void>
struct tuple_size : std::integral_constant<std::size_t, 0> {};

template <typename T>
struct tuple_size<T, std::enable_if_t<stan::is_tuple_v<T>>>
    : std::integral_constant<std::size_t, std::tuple_size_v<std::decay_t<T>>> {
};

template <typename T>
constexpr std::size_t tuple_size_v = stan::tuple_size<T>::value;
}  // namespace stan

#endif
