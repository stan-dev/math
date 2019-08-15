#ifndef STAN_MATH_PRIM_SCAL_META_ENABLE_IF_STAN_SCALAR_HPP
#define STAN_MATH_PRIM_SCAL_META_ENABLE_IF_STAN_SCALAR_HPP

#include <stan/math/prim/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/is_fvar.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>
#include <type_traits>

namespace stan {

template<typename T>
struct is_stan_scalar : std::integral_constant<bool,
 std::is_arithmetic<std::decay_t<T>>::value || is_var<std::decay_t<T>>::value || is_fvar<std::decay_t<T>>::value> {};

template <typename T>
using enable_if_stan_scalar
    = std::enable_if_t<is_stan_scalar<T>::value>;

template <typename T>
using enable_if_not_stan_scalar
    = std::enable_if_t<is_stan_scalar<T>::value>;

template <typename... Types>
using enable_if_all_stan_scalar = std::enable_if_t<
    math::conjunction<is_stan_scalar<std::decay_t<Types>>...>::value>;

template <typename... Types>
using enable_if_any_stan_scalar = std::enable_if_t<
    math::disjunction<is_stan_scalar<std::decay_t<Types>>...>::value>;

template <typename... Types>
using enable_if_all_not_stan_scalar = std::enable_if_t<
    !math::conjunction<is_stan_scalar<std::decay_t<Types>>...>::value>;

template <typename... Types>
using enable_if_any_not_stan_stan_scalar = std::enable_if_t<
    !math::disjunction<is_stan_scalar<std::decay_t<Types>>...>::value>;

}  // namespace stan
#endif
