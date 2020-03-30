#ifndef STAN_MATH_PRIM_META_IS_VAR_OR_ARITHMETIC_HPP
#define STAN_MATH_PRIM_META_IS_VAR_OR_ARITHMETIC_HPP

#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/conjunction.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>

#include <type_traits>

namespace stan {

/** \ingroup type_trait
 * Defines a static member value which is defined to be true (1)
 * if the unqualified cv of type T or its underlying type (if a container)
 * is either var or an arithmetic type, and false (0) otherwise.
 */
template <typename T>
struct is_var_or_arithmetic_type
    : bool_constant<(
          is_var<scalar_type_t<std::decay_t<T>>>::value
          || std::is_arithmetic<scalar_type_t<std::decay_t<T>>>::value)> {};

/** \ingroup type_trait
 * Extends std::true_type if all the provided types are either var or
 * an arithmetic type, extends std::false_type otherwise.
 */
template <typename... T>
using is_var_or_arithmetic = math::conjunction<is_var_or_arithmetic_type<T>...>;

STAN_ADD_REQUIRE_UNARY(var_or_arithmetic, is_var_or_arithmetic,
                       require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY_INNER(var_or_arithmetic, is_var_or_arithmetic,
                             require_stan_scalar_real);

}  // namespace stan
#endif
