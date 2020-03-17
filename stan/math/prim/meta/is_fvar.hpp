#ifndef STAN_MATH_PRIM_META_IS_FVAR_HPP
#define STAN_MATH_PRIM_META_IS_FVAR_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {
/** \ingroup type_trait
 * Defines a static member function type which is defined to be false
 * as the primitive scalar types cannot be a stan::math::fvar type.
 */
template <typename T, typename = void>
struct is_fvar : std::false_type {};

STAN_ADD_REQUIRE_UNARY(fvar, is_fvar);

}  // namespace stan
#endif
