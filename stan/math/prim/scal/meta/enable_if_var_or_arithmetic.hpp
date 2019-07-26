#ifndef STAN_MATH_PRIM_SCAL_META_ENABLE_IF_VAR_OR_ARITHMETIC_HPP
#define STAN_MATH_PRIM_SCAL_META_ENABLE_IF_VAR_OR_ARITHMETIC_HPP

#include <stan/math/prim/scal/meta/is_var_or_arithmetic.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>

#include <type_traits>

namespace stan {

template <typename T>
using enable_if_var_or_arithmetic
    = std::enable_if_t<is_var_or_arithmetic<T>::value>;

template <typename T>
using enable_if_not_var_or_arithmetic
    = std::enable_if_t<!is_var_or_arithmetic<T>::value>;

template <typename... Types>
using enable_if_all_var_or_arithmetic = std::enable_if_t<
    math::conjunction<is_var_or_arithmetic<Types>...>::value>;

template <typename... Types>
using enable_if_any_var_or_arithmetic = std::enable_if_t<
    math::disjunction<is_var_or_arithmetic<Types>...>::value>;

template <typename... Types>
using enable_if_all_not_var_or_arithmetic = std::enable_if_t<
    !math::conjunction<is_var_or_arithmetic<Types>...>::value>;

template <typename... Types>
using enable_if_any_not_var_or_arithmetic = std::enable_if_t<
    !math::disjunction<is_var_or_arithmetic<Types>...>::value>;

}  // namespace stan
#endif
