#ifndef STAN_MATH_PRIM_META_IS_VAR_HPP
#define STAN_MATH_PRIM_META_IS_VAR_HPP

#include <stan/math/prim/meta/require_helpers.hpp>

#include <type_traits>

namespace stan {
/** \ingroup type_trait
 * Defines a static member named value which is defined to be false
 * as the primitive scalar types cannot be a stan::math::var type.
 */
template <typename T, typename = void>
struct is_var : std::false_type {};

/** \addtogroup require_stan_scalar
*  @{
*/
STAN_ADD_REQUIRE_UNARY(var, is_var);
STAN_ADD_REQUIRE_UNARY_SCALAR(var, is_var);
STAN_ADD_REQUIRE_UNARY_VALUE(var, is_var);
/** @}*/
}  // namespace stan
#endif
