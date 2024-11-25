#ifndef STAN_MATH_PRIM_META_IS_ALL_ARITHMETIC_HPP
#define STAN_MATH_PRIM_META_IS_ALL_ARITHMETIC_HPP

#include <stan/math/prim/meta/scalar_type.hpp>
#include <type_traits>
#include <tuple>

namespace stan {
template <typename T>
using has_arithmetic_scalar_type = std::is_arithmetic<scalar_type_t<T>>;

namespace internal {

template <typename... Types>
struct is_all_arithmetic_scalar_impl : std::conjunction<has_arithmetic_scalar_type<std::decay_t<Types>>...> {};

template <typename... Types>
struct is_all_arithmetic_scalar_impl<std::tuple<Types...>> : std::conjunction<is_all_arithmetic_scalar_impl<scalar_type_t<std::decay_t<Types>>>...> {};
}

template <typename... Types>
struct is_all_arithmetic_scalar : std::conjunction<internal::is_all_arithmetic_scalar_impl<std::decay_t<Types>>...> {};

}

#endif