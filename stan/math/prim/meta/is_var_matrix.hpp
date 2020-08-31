#ifndef STAN_MATH_PRIM_META_IS_VAR_MATRIX_HPP
#define STAN_MATH_PRIM_META_IS_VAR_MATRIX_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {
/** \ingroup type_trait
 * Defines a static member named value which is defined to be true
 * if the type is either a type derived from `Eigen::EigenBase` with a `Scalar`
 *  type of `var_value<double>` or a `var_value<T>` where T is derived from
 * `Eigen::EigenBase`
 */
template <typename T, typename = void>
struct is_var_matrix : std::false_type {};

STAN_ADD_REQUIRE_UNARY(var_matrix, is_var_matrix, require_eigens_types);
STAN_ADD_REQUIRE_CONTAINER(var_matrix, is_var_matrix, require_eigens_types);
STAN_ADD_REQUIRE_UNARY_INNER(var_matrix, is_var_matrix, require_eigens_types);

}  // namespace stan
#endif
