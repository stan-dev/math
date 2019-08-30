#ifndef STAN_MATH_PRIM_SCAL_META_ENABLE_IF_DOUBLE_OR_INT_HPP
#define STAN_MATH_PRIM_SCAL_META_ENABLE_IF_DOUBLE_OR_INT_HPP

#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>

#include <type_traits>

namespace stan {

template <typename T>
struct is_double_or_int : std::integral_constant<bool, std::is_same<double, std::decay_t<T>>::value || std::is_same<int, std::decay_t<T>>::value> {};

template <typename T>
using enable_if_double_or_int
    = std::enable_if_t<is_double_or_int<T>::value>;

template <typename T>
using enable_if_not_double_or_int
    = std::enable_if_t<!is_double_or_int<T>::value>;

template <typename... Types>
using enable_if_all_double_or_int = std::enable_if_t<
    math::conjunction<is_double_or_int<Types>...>::value>;

template <typename... Types>
using enable_if_any_double_or_int = std::enable_if_t<
    math::disjunction<is_double_or_int<Types>...>::value>;

template <typename... Types>
using enable_if_all_not_double_or_int = std::enable_if_t<
    !math::conjunction<is_double_or_int<Types>...>::value>;

template <typename... Types>
using enable_if_any_not_double_or_int = std::enable_if_t<
    !math::disjunction<is_double_or_int<Types>...>::value>;

template <typename T>
using double_or_int_type = enable_if_double_or_int<std::decay_t<T>>;

template <typename T>
using not_double_or_int_type = enable_if_not_double_or_int<std::decay_t<T>>;

template <typename... Types>
using all_double_or_int_type = enable_if_all_double_or_int<std::decay_t<Types>...>;

template <typename... Types>
using any_double_or_int_type = enable_if_any_double_or_int<std::decay_t<Types>...>;

template <typename... Types>
using not_all_double_or_int_type = enable_if_all_not_double_or_int<std::decay_t<Types>...>;

template <typename... Types>
using not_any_double_or_int_type = enable_if_any_not_double_or_int<std::decay_t<Types>...>;

}  // namespace stan
#endif
