#ifndef STAN_MATH_PRIM_SCAL_META_ENABLE_IF_contains_var_HPP
#define STAN_MATH_PRIM_SCAL_META_ENABLE_IF_contains_var_HPP

#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <stan/math/prim/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>
#include <type_traits>

namespace stan {

template<typename T>
struct is_contains_var: std::integral_constant<bool, is_var<scalar_type_t<std::decay_t<T>>>::value> {};

template <typename T>
using enable_if_contains_var
    = std::enable_if_t<is_contains_var<T>::value>;

template <typename T>
using enable_if_not_contains_var
    = std::enable_if_t<is_contains_var<T>::value>;

template <typename... Types>
using enable_if_all_contains_var= std::enable_if_t<
    math::conjunction<is_contains_var<Types>...>::value>;

template <typename... Types>
using enable_if_any_contains_var= std::enable_if_t<
    math::disjunction<is_contains_var<Types>...>::value>;

template <typename... Types>
using enable_if_all_not_contains_var= std::enable_if_t<
    !math::conjunction<is_contains_var<Types>...>::value>;

template <typename... Types>
using enable_if_any_not_contains_var= std::enable_if_t<
    !math::disjunction<is_contains_var<Types>...>::value>;

}  // namespace stan
#endif
