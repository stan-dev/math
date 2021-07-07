#ifndef STAN_MATH_PRIM_META_IS_FVAR_OR_ARITHMETIC_HPP
#define STAN_MATH_PRIM_META_IS_FVAR_OR_ARITHMETIC_HPP

#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_fvar.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/conjunction.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>

#include <type_traits>

namespace stan {

/** \ingroup type_trait
 * Extends std::true_type if all the provided types are either fvar or
 * an arithmetic type, extends std::false_type otherwise.
 */
template <typename T>
using is_fvar_or_arithmetic
    = bool_constant<std::is_arithmetic<std::decay_t<T>>::value
                    || is_fvar<std::decay_t<T>>::value>;

STAN_ADD_REQUIRE_UNARY(fvar_or_arithmetic, is_fvar_or_arithmetic,
                       require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY_INNER(fvar_or_arithmetic, is_fvar_or_arithmetic,
                             require_stan_scalar_real);

}  // namespace stan
#endif
