#ifndef STAN_MATH_PRIM_META_IS_MATRIX_HPP
#define STAN_MATH_PRIM_META_IS_MATRIX_HPP

#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/prim/meta/is_rev_matrix.hpp>

namespace stan {
/**
 * Check if a type is derived from `Eigen::EigenBase` or is a `var_value`
 *  whose `value_type` is derived from `Eigen::EigenBase`
 * @tparam T type to check.
 * @ingroup type_trait
 */
template <typename T>
struct is_matrix
    : bool_constant<math::disjunction<is_rev_matrix<T>, is_eigen<T>>::value> {};
STAN_ADD_REQUIRE_UNARY(matrix, is_matrix, require_eigens_types);
STAN_ADD_REQUIRE_UNARY_INNER(matrix, is_matrix, require_eigens_types);
STAN_ADD_REQUIRE_CONTAINER(matrix, is_matrix, require_eigens_types);

}  // namespace stan

#endif
