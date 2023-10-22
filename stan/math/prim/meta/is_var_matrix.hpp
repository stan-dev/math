#ifndef STAN_MATH_PRIM_META_IS_VAR_MATRIX
#define STAN_MATH_PRIM_META_IS_VAR_MATRIX

#include <stan/math/prim/meta/conjunction.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_vector.hpp>

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

//STAN_ADD_REQUIRE_UNARY(var_matrix, is_var_matrix, require_eigens_types);
template <typename T>
using require_var_matrix_t = require_t<is_var_matrix<std::decay_t<T>>>;

template <typename T>
using require_not_var_matrix_t
    = require_not_t<is_var_matrix<std::decay_t<T>>>;

template <typename... Types>
using require_all_var_matrix_t
    = require_all_t<is_var_matrix<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_var_matrix_t
    = require_any_t<is_var_matrix<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_var_matrix_t
    = require_all_not_t<is_var_matrix<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_var_matrix_t
    = require_any_not_t<is_var_matrix<std::decay_t<Types>>...>;

  
//STAN_ADD_REQUIRE_UNARY_INNER(var_matrix, is_var_matrix, require_eigens_types);
template <typename T>
using require_vt_var_matrix
    = require_t<is_var_matrix<value_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_vt_var_matrix
    = require_not_t<is_var_matrix<value_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_vt_var_matrix
    = require_all_t<is_var_matrix<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_vt_var_matrix
    = require_any_t<is_var_matrix<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_vt_var_matrix
    = require_all_not_t<is_var_matrix<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_vt_var_matrix
    = require_any_not_t<is_var_matrix<value_type_t<std::decay_t<Types>>>...>;

template <typename T>
using require_st_var_matrix
    = require_t<is_var_matrix<scalar_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_st_var_matrix
    = require_not_t<is_var_matrix<scalar_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_st_var_matrix
    = require_all_t<is_var_matrix<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_st_var_matrix
    = require_any_t<is_var_matrix<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_st_var_matrix
    = require_all_not_t<is_var_matrix<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_st_var_matrix
    = require_any_not_t<is_var_matrix<scalar_type_t<std::decay_t<Types>>>...>;


/**
 * Check if a type is a `var_value` whose `value_type` is derived from
 * `Eigen::EigenBase`. And the type must have a compile time constant number
 *  of columns equal to 1.
 * @tparam T type to check.
 * @ingroup type_trait
 */
template <typename T>
struct is_var_col_vector
    : bool_constant<math::conjunction<
          is_var<T>, is_eigen_col_vector<value_type_t<T>>>::value> {};

//STAN_ADD_REQUIRE_UNARY(var_col_vector, is_var_col_vector, require_eigens_types);
template <typename T>
using require_var_col_vector_t = require_t<is_var_col_vector<std::decay_t<T>>>;

template <typename T>
using require_not_var_col_vector_t
    = require_not_t<is_var_col_vector<std::decay_t<T>>>;

template <typename... Types>
using require_all_var_col_vector_t
    = require_all_t<is_var_col_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_var_col_vector_t
    = require_any_t<is_var_col_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_var_col_vector_t
    = require_all_not_t<is_var_col_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_var_col_vector_t
    = require_any_not_t<is_var_col_vector<std::decay_t<Types>>...>;

  
//STAN_ADD_REQUIRE_UNARY_INNER(var_col_vector, is_var_col_vector, require_eigens_types);
template <typename T>
using require_vt_var_col_vector
    = require_t<is_var_col_vector<value_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_vt_var_col_vector
    = require_not_t<is_var_col_vector<value_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_vt_var_col_vector
    = require_all_t<is_var_col_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_vt_var_col_vector
    = require_any_t<is_var_col_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_vt_var_col_vector
    = require_all_not_t<is_var_col_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_vt_var_col_vector
    = require_any_not_t<is_var_col_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename T>
using require_st_var_col_vector
    = require_t<is_var_col_vector<scalar_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_st_var_col_vector
    = require_not_t<is_var_col_vector<scalar_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_st_var_col_vector
    = require_all_t<is_var_col_vector<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_st_var_col_vector
    = require_any_t<is_var_col_vector<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_st_var_col_vector
    = require_all_not_t<is_var_col_vector<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_st_var_col_vector
    = require_any_not_t<is_var_col_vector<scalar_type_t<std::decay_t<Types>>>...>;


/**
 * Check if a type is a `var_value` whose `value_type` is derived from
 * `Eigen::EigenBase`. And the type must have a compile time constant number
 *  of rows equal to 1.
 * @tparam T type to check.
 * @ingroup type_trait
 */
template <typename T>
struct is_var_row_vector
    : bool_constant<math::conjunction<
          is_var<T>, is_eigen_row_vector<value_type_t<T>>>::value> {};

//STAN_ADD_REQUIRE_UNARY(var_row_vector, is_var_row_vector, require_eigens_types);
template <typename T>
using require_var_row_vector_t = require_t<is_var_row_vector<std::decay_t<T>>>;

template <typename T>
using require_not_var_row_vector_t
    = require_not_t<is_var_row_vector<std::decay_t<T>>>;

template <typename... Types>
using require_all_var_row_vector_t
    = require_all_t<is_var_row_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_var_row_vector_t
    = require_any_t<is_var_row_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_var_row_vector_t
    = require_all_not_t<is_var_row_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_var_row_vector_t
    = require_any_not_t<is_var_row_vector<std::decay_t<Types>>...>;

  
//STAN_ADD_REQUIRE_UNARY_INNER(var_row_vector, is_var_row_vector, require_eigens_types);
template <typename T>
using require_vt_var_row_vector
    = require_t<is_var_row_vector<value_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_vt_var_row_vector
    = require_not_t<is_var_row_vector<value_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_vt_var_row_vector
    = require_all_t<is_var_row_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_vt_var_row_vector
    = require_any_t<is_var_row_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_vt_var_row_vector
    = require_all_not_t<is_var_row_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_vt_var_row_vector
    = require_any_not_t<is_var_row_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename T>
using require_st_var_row_vector
    = require_t<is_var_row_vector<scalar_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_st_var_row_vector
    = require_not_t<is_var_row_vector<scalar_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_st_var_row_vector
    = require_all_t<is_var_row_vector<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_st_var_row_vector
    = require_any_t<is_var_row_vector<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_st_var_row_vector
    = require_all_not_t<is_var_row_vector<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_st_var_row_vector
    = require_any_not_t<is_var_row_vector<scalar_type_t<std::decay_t<Types>>>...>;


/**
 * Check if a type is a `var_value` whose `value_type` is derived from
 * `Eigen::EigenBase`. And the type must have a compile time constant number
 *  of columns or rows equal to 1.
 * @tparam T type to check.
 * @ingroup type_trait
 */
template <typename T>
struct is_var_vector
    : bool_constant<math::disjunction<is_var_col_vector<T>,
                                      is_var_row_vector<T>>::value> {};

//STAN_ADD_REQUIRE_UNARY(var_vector, is_var_vector, require_eigens_types);
template <typename T>
using require_var_vector_t = require_t<is_var_vector<std::decay_t<T>>>;

template <typename T>
using require_not_var_vector_t
    = require_not_t<is_var_vector<std::decay_t<T>>>;

template <typename... Types>
using require_all_var_vector_t
    = require_all_t<is_var_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_var_vector_t
    = require_any_t<is_var_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_var_vector_t
    = require_all_not_t<is_var_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_var_vector_t
    = require_any_not_t<is_var_vector<std::decay_t<Types>>...>;

  
//STAN_ADD_REQUIRE_UNARY_INNER(var_vector, is_var_vector, require_eigens_types);
template <typename T>
using require_vt_var_vector
    = require_t<is_var_vector<value_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_vt_var_vector
    = require_not_t<is_var_vector<value_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_vt_var_vector
    = require_all_t<is_var_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_vt_var_vector
    = require_any_t<is_var_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_vt_var_vector
    = require_all_not_t<is_var_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_vt_var_vector
    = require_any_not_t<is_var_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename T>
using require_st_var_vector
    = require_t<is_var_vector<scalar_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_st_var_vector
    = require_not_t<is_var_vector<scalar_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_st_var_vector
    = require_all_t<is_var_vector<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_st_var_vector
    = require_any_t<is_var_vector<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_st_var_vector
    = require_all_not_t<is_var_vector<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_st_var_vector
    = require_any_not_t<is_var_vector<scalar_type_t<std::decay_t<Types>>>...>;


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
