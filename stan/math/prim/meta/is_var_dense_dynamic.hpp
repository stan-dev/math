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

//STAN_ADD_REQUIRE_UNARY(var_dense_dynamic, is_var_dense_dynamic, require_eigens_types);
template <typename T>
using require_var_dense_dynamic_t = require_t<is_var_dense_dynamic<std::decay_t<T>>>;

template <typename T>
using require_not_var_dense_dynamic_t
    = require_not_t<is_var_dense_dynamic<std::decay_t<T>>>;

template <typename... Types>
using require_all_var_dense_dynamic_t
    = require_all_t<is_var_dense_dynamic<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_var_dense_dynamic_t
    = require_any_t<is_var_dense_dynamic<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_var_dense_dynamic_t
    = require_all_not_t<is_var_dense_dynamic<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_var_dense_dynamic_t
    = require_any_not_t<is_var_dense_dynamic<std::decay_t<Types>>...>;

  
//STAN_ADD_REQUIRE_UNARY_INNER(var_dense_dynamic, is_var_dense_dynamic, require_eigens_types);
template <typename T>
using require_vt_var_dense_dynamic
    = require_t<is_var_dense_dynamic<value_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_vt_var_dense_dynamic
    = require_not_t<is_var_dense_dynamic<value_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_vt_var_dense_dynamic
    = require_all_t<is_var_dense_dynamic<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_vt_var_dense_dynamic
    = require_any_t<is_var_dense_dynamic<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_vt_var_dense_dynamic
    = require_all_not_t<is_var_dense_dynamic<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_vt_var_dense_dynamic
    = require_any_not_t<is_var_dense_dynamic<value_type_t<std::decay_t<Types>>>...>;

template <typename T>
using require_st_var_dense_dynamic
    = require_t<is_var_dense_dynamic<scalar_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_st_var_dense_dynamic
    = require_not_t<is_var_dense_dynamic<scalar_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_st_var_dense_dynamic
    = require_all_t<is_var_dense_dynamic<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_st_var_dense_dynamic
    = require_any_t<is_var_dense_dynamic<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_st_var_dense_dynamic
    = require_all_not_t<is_var_dense_dynamic<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_st_var_dense_dynamic
    = require_any_not_t<is_var_dense_dynamic<scalar_type_t<std::decay_t<Types>>>...>;


}  // namespace stan

#endif
