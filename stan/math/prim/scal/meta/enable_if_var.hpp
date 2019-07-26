#ifndef STAN_MATH_PRIM_SCAL_META_ENABLE_IF_VAR_HPP
#define STAN_MATH_PRIM_SCAL_META_ENABLE_IF_VAR_HPP

#include <stan/math/prim/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>

#include <type_traits>

namespace stan {

template <typename T>
using enable_if_var = std::enable_if_t<is_var<T>::value>;

template <typename... Types>
using enable_if_all_var
    = std::enable_if_t<math::conjunction<is_var<Types>...>::value>;

template <typename... Types>
using enable_if_any_var
    = std::enable_if_t<math::disjunction<is_var<Types>...>::value>;

template <typename T>
using enable_if_not_var = std::enable_if_t<!is_var<T>::value>;

template <typename... Types>
using enable_if_all_not_var
    = std::enable_if_t<!math::conjunction<is_var<Types>...>::value>;

template <typename... Types>
using enable_if_any_not_var
    = std::enable_if_t<!math::disjunction<is_var<Types>...>::value>;

}  // namespace stan
#endif
