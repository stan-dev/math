#ifndef STAN_MATH_PRIM_SCAL_META_ENABLE_IF_SCALAR_HPP
#define STAN_MATH_PRIM_SCAL_META_ENABLE_IF_SCALAR_HPP

#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>

#include <type_traits>

namespace stan {

template <typename T>
using enable_if_scalar = std::enable_if_t<std::is_scalar<T>::value>;

template <typename T>
using enable_if_not_scalar = std::enable_if_t<!std::is_scalar<T>::value>;

template <typename... Types>
using enable_if_all_scalar
    = std::enable_if_t<math::conjunction<std::is_scalar<Types>...>::value>;

template <typename... Types>
using enable_if_any_scalar
    = std::enable_if_t<math::disjunction<std::is_scalar<Types>...>::value>;

template <typename... Types>
using enable_if_all_not_scalar
    = std::enable_if_t<!math::conjunction<std::is_scalar<Types>...>::value>;

template <typename... Types>
using enable_if_any_not_scalar
    = std::enable_if_t<!math::disjunction<std::is_scalar<Types>...>::value>;

}  // namespace stan
#endif
