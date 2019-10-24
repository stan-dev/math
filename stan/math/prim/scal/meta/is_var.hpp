#ifndef STAN_MATH_PRIM_SCAL_META_IS_VAR_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_VAR_HPP

#include <type_traits>

namespace stan {
/**
 * Defines a static member named value which is defined to be false
 * as the primitive scalar types cannot be a stan::math::var type.
 */
template <typename T, typename = void>
struct is_var : std::false_type {};

}  // namespace stan
#endif
