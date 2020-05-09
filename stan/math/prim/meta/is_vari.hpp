#ifndef STAN_MATH_PRIM_META_IS_VARI_HPP
#define STAN_MATH_PRIM_META_IS_VARI_HPP

#include <stan/math/prim/meta/require_helpers.hpp>

#include <type_traits>

namespace stan {
/** \ingroup type_trait
 * Defines a static member named value which is defined to be false
 * as the primitive scalar types cannot be a stan::math::var type.
 */
template <typename T, typename = void>
struct is_vari : std::false_type {};

STAN_ADD_REQUIRE_UNARY(vari, is_vari, require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY_INNER(vari, is_vari, require_stan_scalar_real);

}  // namespace stan
#endif
