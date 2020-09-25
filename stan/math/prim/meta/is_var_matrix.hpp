#ifndef STAN_MATH_PRIM_META_IS_VAR_MATRIX
#define STAN_MATH_PRIM_META_IS_VAR_MATRIX

#include <stan/math/prim/meta/conjunction.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>

namespace stan {
/**
 * Check if a type is a `var_value` whose `value_type` is derived from
 * `Eigen::EigenBase`
 * @tparam T type to check.
 * @ingroup type_trait
 */
template <typename T>
struct is_var_matrix
    : bool_constant<
          math::conjunction<is_var<T>, is_eigen<value_type_t<T>>>::value> {};

STAN_ADD_REQUIRE_UNARY(var_matrix, is_var_matrix, require_eigens_types);
STAN_ADD_REQUIRE_UNARY_INNER(var_matrix, is_var_matrix, require_eigens_types);

/**
 * Check if any types in a parameter pack are a `var_value` whose `value_type`
 *  is derived from `Eigen::EigenBase`
 * @tparam Types parameter pack of types to check.
 * @ingroup type_trait
 */
template <typename... Types>
struct is_any_var_matrix
    : bool_constant<math::disjunction<is_var_matrix<Types>...>::value> {};

}  // namespace stan

#endif
