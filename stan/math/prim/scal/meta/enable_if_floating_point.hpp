#ifndef STAN_MATH_PRIM_SCAL_META_ENABLE_IF_FLOATING_POINT_HPP
#define STAN_MATH_PRIM_SCAL_META_ENABLE_IF_FLOATING_POINT_HPP

#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>

#include <type_traits>

namespace stan {

template <typename T>
using enable_if_floating_point
    = std::enable_if_t<std::is_floating_point<T>::value>;

template <typename T>
using enable_if_not_floating_point
    = std::enable_if_t<!std::is_floating_point<T>::value>;

template <typename... Types>
using enable_if_all_floating_point = std::enable_if_t<
    math::conjunction<std::is_floating_point<Types>...>::value>;

template <typename... Types>
using enable_if_any_floating_point = std::enable_if_t<
    math::disjunction<std::is_floating_point<Types>...>::value>;

template <typename... Types>
using enable_if_all_not_floating_point = std::enable_if_t<
    !math::conjunction<std::is_floating_point<Types>...>::value>;

template <typename... Types>
using enable_if_any_not_floating_point = std::enable_if_t<
    !math::disjunction<std::is_floating_point<Types>...>::value>;

}  // namespace stan
#endif
