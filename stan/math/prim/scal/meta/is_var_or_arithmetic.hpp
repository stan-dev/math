#ifndef STAN_MATH_PRIM_SCAL_META_IS_VAR_OR_ARITHMETIC_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_VAR_OR_ARITHMETIC_HPP

#include <stan/math/prim/scal/meta/bool_constant.hpp>
#include <stan/math/prim/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <type_traits>

namespace stan {

/**
 * Defines a static member value which is defined to be true (1)
 * if the type is either var or an aritmetic type
 * and false (0) otherwise.
 */
template <typename T>
struct is_var_or_arithmetic_type
    : bool_constant<(is_var<scalar_type_t<T>>::value
                     || std::is_arithmetic<scalar_type_t<T>>::value)> {};

/**
 * Extends std::true_type if all the provided types are either var or
 * an arithmetic type, extends std::false_type otherwise.
 */
template <typename... T>
using is_var_or_arithmetic = math::conjunction<is_var_or_arithmetic_type<T>...>;

template<class... T>
constexpr bool is_var_or_arithmetic_v = is_var_or_arithmetic<T...>::value;

}  // namespace stan
#endif
