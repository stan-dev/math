#ifndef STAN_MATH_PRIM_SCAL_META_ENABLE_IF_CONTAINS_STAN_SCALAR_HPP
#define STAN_MATH_PRIM_SCAL_META_ENABLE_IF_CONTAINS_STAN_SCALAR_HPP

#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <stan/math/prim/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/is_fvar.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>
#include <type_traits>

namespace stan {

template <typename T>
struct is_contains_stan_scalar
    : std::integral_constant<
          bool, std::is_arithmetic<scalar_type_t<std::decay_t<T>>>::value
                    || is_var<scalar_type_t<std::decay_t<T>>>::value
                    || is_fvar<scalar_type_t<std::decay_t<T>>>::value> {};

template <typename T>
using enable_if_contains_stan_scalar
    = std::enable_if_t<is_contains_stan_scalar<T>::value>;

template <typename T>
using enable_if_not_contains_stan_scalar
    = std::enable_if_t<is_contains_stan_scalar<T>::value>;

template <typename... Types>
using enable_if_all_contains_stan_scalar = std::enable_if_t<
    math::conjunction<is_contains_stan_scalar<Types>...>::value>;

template <typename... Types>
using enable_if_any_contains_stan_scalar = std::enable_if_t<
    math::disjunction<is_contains_stan_scalar<Types>...>::value>;

template <typename... Types>
using enable_if_all_not_contains_stan_scalar = std::enable_if_t<
    !math::conjunction<is_contains_stan_scalar<Types>...>::value>;

template <typename... Types>
using enable_if_any_not_stan_contains_stan_scalar = std::enable_if_t<
    !math::disjunction<is_contains_stan_scalar<Types>...>::value>;

}  // namespace stan
#endif
