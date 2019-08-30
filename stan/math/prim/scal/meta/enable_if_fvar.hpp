#ifndef STAN_MATH_PRIM_SCAL_META_ENABLE_IF_FVAR_HPP
#define STAN_MATH_PRIM_SCAL_META_ENABLE_IF_FVAR_HPP

#include <stan/math/prim/scal/meta/is_fvar.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>

#include <type_traits>

namespace stan {

template <typename T>
using enable_if_fvar = std::enable_if_t<is_fvar<T>::value>;

template <typename... Types>
using enable_if_all_fvar
    = std::enable_if_t<math::conjunction<is_fvar<Types>...>::value>;

template <typename... Types>
using enable_if_any_fvar
    = std::enable_if_t<math::disjunction<is_fvar<Types>...>::value>;

template <typename T>
using enable_if_not_fvar = std::enable_if_t<!is_fvar<T>::value>;

template <typename... Types>
using enable_if_all_not_fvar
    = std::enable_if_t<!math::conjunction<is_fvar<Types>...>::value>;

template <typename... Types>
using enable_if_any_not_fvar
    = std::enable_if_t<!math::disjunction<is_fvar<Types>...>::value>;

}  // namespace stan
#endif
