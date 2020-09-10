#ifndef STAN_MATH_PRIM_META_IS_PLAIN_TYPE_HPP
#define STAN_MATH_PRIM_META_IS_PLAIN_TYPE_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <stan/math/prim/meta/plain_type.hpp>
#include <type_traits>

namespace stan {
/** \ingroup type_trait
 * Checks whether the template type `T` is an assignable type. This is used
 * to detect whether a type is an Eigen matrix expression.
 */
template <typename S>
using is_plain_type = std::is_same<std::decay_t<S>, plain_type_t<S>>;

STAN_ADD_REQUIRE_UNARY(plain_type, is_plain_type, require_eigens_types);
STAN_ADD_REQUIRE_UNARY_INNER(plain_type, is_plain_type, require_eigens_types);

}  // namespace stan
#endif
