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

STAN_ADD_REQUIRE_UNARY(var, is_var, require_stan_scalar_real);
STAN_ADD_REQUIRE_CONTAINER(var, is_var, require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY_INNER(var, is_var, require_stan_scalar_real);

template <typename T>
struct value_type<T, std::enable_if_t<is_var<T>::value>> {
  using type = typename std::decay_t<T>::value_type;
};

}  // namespace stan
#endif
