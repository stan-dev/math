#ifndef STAN_MATH_PRIM_META_IS_VAR_DENSE_DYNAMIC
#define STAN_MATH_PRIM_META_IS_VAR_DENSE_DYNAMIC

#include <stan/math/prim/meta/conjunction.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/prim/meta/is_eigen_dense_dynamic.hpp>

namespace stan {
/**
 * Check if a type is a `var_value` whose `value_type` is derived from
 * `Eigen::EigenBase` and has dynamic rows and columns
 * @tparam T type to check.
 * @ingroup type_trait
 */
template <typename T>
struct is_var_dense_dynamic
    : bool_constant<math::conjunction<
          is_var<T>, is_eigen_dense_dynamic<value_type_t<T>>>::value> {};

STAN_ADD_REQUIRE_UNARY(var_dense_dynamic, is_var_dense_dynamic,
                       require_eigens_types);
STAN_ADD_REQUIRE_UNARY_INNER(var_dense_dynamic, is_var_dense_dynamic,
                             require_eigens_types);

}  // namespace stan

#endif
