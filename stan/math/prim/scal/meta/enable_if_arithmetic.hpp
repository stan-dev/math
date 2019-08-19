#ifndef STAN_MATH_PRIM_SCAL_META_ENABLE_IF_ARITHMETIC_HPP
#define STAN_MATH_PRIM_SCAL_META_ENABLE_IF_ARITHMETIC_HPP

#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>

#include <type_traits>

namespace stan {

template <typename T>
using enable_if_arithmetic = std::enable_if_t<std::is_arithmetic<T>::value>;

template <typename T>
using enable_if_not_arithmetic
    = std::enable_if_t<!std::is_arithmetic<T>::value>;

template <typename... Types>
using enable_if_all_arithmetic
    = std::enable_if_t<math::conjunction<std::is_arithmetic<Types>...>::value>;

template <typename... Types>
using enable_if_any_arithmetic
    = std::enable_if_t<math::disjunction<std::is_arithmetic<Types>...>::value>;

template <typename... Types>
using enable_if_all_not_arithmetic
    = std::enable_if_t<!math::conjunction<std::is_arithmetic<Types>...>::value>;

template <typename... Types>
using enable_if_any_not_arithmetic
    = std::enable_if_t<!math::disjunction<std::is_arithmetic<Types>...>::value>;

}  // namespace stan
#endif
