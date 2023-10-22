#ifndef STAN_MATH_PRIM_META_IS_REV_MATRIX_HPP
#define STAN_MATH_PRIM_META_IS_REV_MATRIX_HPP

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
struct is_rev_matrix : std::false_type {};

//STAN_ADD_REQUIRE_UNARY(rev_matrix, is_rev_matrix, require_eigens_types);
template <typename T>
using require_rev_matrix_t = require_t<is_rev_matrix<std::decay_t<T>>>;

template <typename T>
using require_not_rev_matrix_t
    = require_not_t<is_rev_matrix<std::decay_t<T>>>;

template <typename... Types>
using require_all_rev_matrix_t
    = require_all_t<is_rev_matrix<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_rev_matrix_t
    = require_any_t<is_rev_matrix<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_rev_matrix_t
    = require_all_not_t<is_rev_matrix<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_rev_matrix_t
    = require_any_not_t<is_rev_matrix<std::decay_t<Types>>...>;

  
//STAN_ADD_REQUIRE_CONTAINER(rev_matrix, is_rev_matrix, require_eigens_types);
template <template <class...> class TypeCheck, class... Check>
using require_rev_matrix_vt = require_t<
    container_type_check_base<is_rev_matrix, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_rev_matrix_vt = require_not_t<
    container_type_check_base<is_rev_matrix, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_rev_matrix_vt = require_any_t<
    container_type_check_base<is_rev_matrix, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_rev_matrix_vt = require_any_not_t<
    container_type_check_base<is_rev_matrix, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_rev_matrix_vt = require_all_t<
    container_type_check_base<is_rev_matrix, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_rev_matrix_vt = require_all_not_t<
    container_type_check_base<is_rev_matrix, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_rev_matrix_st = require_t<
    container_type_check_base<is_rev_matrix, scalar_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_rev_matrix_st = require_not_t<
    container_type_check_base<is_rev_matrix, scalar_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_rev_matrix_st = require_any_t<
    container_type_check_base<is_rev_matrix, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_rev_matrix_st = require_any_not_t<
    container_type_check_base<is_rev_matrix, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_rev_matrix_st = require_all_t<
    container_type_check_base<is_rev_matrix, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_rev_matrix_st = require_all_not_t<
    container_type_check_base<is_rev_matrix, scalar_type_t, TypeCheck, Check>...>;

//STAN_ADD_REQUIRE_UNARY_INNER(rev_matrix, is_rev_matrix, require_eigens_types);
template <typename T>
using require_vt_rev_matrix
    = require_t<is_rev_matrix<value_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_vt_rev_matrix
    = require_not_t<is_rev_matrix<value_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_vt_rev_matrix
    = require_all_t<is_rev_matrix<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_vt_rev_matrix
    = require_any_t<is_rev_matrix<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_vt_rev_matrix
    = require_all_not_t<is_rev_matrix<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_vt_rev_matrix
    = require_any_not_t<is_rev_matrix<value_type_t<std::decay_t<Types>>>...>;

template <typename T>
using require_st_rev_matrix
    = require_t<is_rev_matrix<scalar_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_st_rev_matrix
    = require_not_t<is_rev_matrix<scalar_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_st_rev_matrix
    = require_all_t<is_rev_matrix<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_st_rev_matrix
    = require_any_t<is_rev_matrix<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_st_rev_matrix
    = require_all_not_t<is_rev_matrix<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_st_rev_matrix
    = require_any_not_t<is_rev_matrix<scalar_type_t<std::decay_t<Types>>>...>;


/** \ingroup type_trait
 * Defines a static member named value which is defined to be true
 * if the type is either derived from `Eigen::EigenBase` with a `Scalar`
 *  type of `var_value<double>` or a `var_value<T>` where T is derived from
 * `Eigen::EigenBase`. And the type must have a compile time constant number
 *  of columns equal to 1.
 */
template <typename T, typename = void>
struct is_rev_col_vector : std::false_type {};

//STAN_ADD_REQUIRE_UNARY(rev_col_vector, is_rev_col_vector, require_eigens_types);
template <typename T>
using require_rev_col_vector_t = require_t<is_rev_col_vector<std::decay_t<T>>>;

template <typename T>
using require_not_rev_col_vector_t
    = require_not_t<is_rev_col_vector<std::decay_t<T>>>;

template <typename... Types>
using require_all_rev_col_vector_t
    = require_all_t<is_rev_col_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_rev_col_vector_t
    = require_any_t<is_rev_col_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_rev_col_vector_t
    = require_all_not_t<is_rev_col_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_rev_col_vector_t
    = require_any_not_t<is_rev_col_vector<std::decay_t<Types>>...>;

  
//STAN_ADD_REQUIRE_CONTAINER(rev_col_vector, is_rev_col_vector, require_eigens_types);
template <template <class...> class TypeCheck, class... Check>
using require_rev_col_vector_vt = require_t<
    container_type_check_base<is_rev_col_vector, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_rev_col_vector_vt = require_not_t<
    container_type_check_base<is_rev_col_vector, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_rev_col_vector_vt = require_any_t<
    container_type_check_base<is_rev_col_vector, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_rev_col_vector_vt = require_any_not_t<
    container_type_check_base<is_rev_col_vector, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_rev_col_vector_vt = require_all_t<
    container_type_check_base<is_rev_col_vector, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_rev_col_vector_vt = require_all_not_t<
    container_type_check_base<is_rev_col_vector, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_rev_col_vector_st = require_t<
    container_type_check_base<is_rev_col_vector, scalar_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_rev_col_vector_st = require_not_t<
    container_type_check_base<is_rev_col_vector, scalar_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_rev_col_vector_st = require_any_t<
    container_type_check_base<is_rev_col_vector, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_rev_col_vector_st = require_any_not_t<
    container_type_check_base<is_rev_col_vector, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_rev_col_vector_st = require_all_t<
    container_type_check_base<is_rev_col_vector, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_rev_col_vector_st = require_all_not_t<
    container_type_check_base<is_rev_col_vector, scalar_type_t, TypeCheck, Check>...>;

//STAN_ADD_REQUIRE_UNARY_INNER(rev_col_vector, is_rev_col_vector, require_eigens_types);
template <typename T>
using require_vt_rev_col_vector
    = require_t<is_rev_col_vector<value_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_vt_rev_col_vector
    = require_not_t<is_rev_col_vector<value_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_vt_rev_col_vector
    = require_all_t<is_rev_col_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_vt_rev_col_vector
    = require_any_t<is_rev_col_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_vt_rev_col_vector
    = require_all_not_t<is_rev_col_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_vt_rev_col_vector
    = require_any_not_t<is_rev_col_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename T>
using require_st_rev_col_vector
    = require_t<is_rev_col_vector<scalar_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_st_rev_col_vector
    = require_not_t<is_rev_col_vector<scalar_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_st_rev_col_vector
    = require_all_t<is_rev_col_vector<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_st_rev_col_vector
    = require_any_t<is_rev_col_vector<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_st_rev_col_vector
    = require_all_not_t<is_rev_col_vector<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_st_rev_col_vector
    = require_any_not_t<is_rev_col_vector<scalar_type_t<std::decay_t<Types>>>...>;


/** \ingroup type_trait
 * Defines a static member named value which is defined to be true
 * if the type is either a type derived from `Eigen::EigenBase` with a `Scalar`
 *  type of `var_value<double>` or a `var_value<T>` where T is derived from
 * `Eigen::EigenBase`. And the type must have a compile time constant number
 *  of rows equal to 1.
 */
template <typename T, typename = void>
struct is_rev_row_vector : std::false_type {};

//STAN_ADD_REQUIRE_UNARY(rev_row_vector, is_rev_row_vector, require_eigens_types);
template <typename T>
using require_rev_row_vector_t = require_t<is_rev_row_vector<std::decay_t<T>>>;

template <typename T>
using require_not_rev_row_vector_t
    = require_not_t<is_rev_row_vector<std::decay_t<T>>>;

template <typename... Types>
using require_all_rev_row_vector_t
    = require_all_t<is_rev_row_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_rev_row_vector_t
    = require_any_t<is_rev_row_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_rev_row_vector_t
    = require_all_not_t<is_rev_row_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_rev_row_vector_t
    = require_any_not_t<is_rev_row_vector<std::decay_t<Types>>...>;

  
//STAN_ADD_REQUIRE_CONTAINER(rev_row_vector, is_rev_row_vector, require_eigens_types);
template <template <class...> class TypeCheck, class... Check>
using require_rev_row_vector_vt = require_t<
    container_type_check_base<is_rev_row_vector, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_rev_row_vector_vt = require_not_t<
    container_type_check_base<is_rev_row_vector, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_rev_row_vector_vt = require_any_t<
    container_type_check_base<is_rev_row_vector, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_rev_row_vector_vt = require_any_not_t<
    container_type_check_base<is_rev_row_vector, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_rev_row_vector_vt = require_all_t<
    container_type_check_base<is_rev_row_vector, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_rev_row_vector_vt = require_all_not_t<
    container_type_check_base<is_rev_row_vector, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_rev_row_vector_st = require_t<
    container_type_check_base<is_rev_row_vector, scalar_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_rev_row_vector_st = require_not_t<
    container_type_check_base<is_rev_row_vector, scalar_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_rev_row_vector_st = require_any_t<
    container_type_check_base<is_rev_row_vector, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_rev_row_vector_st = require_any_not_t<
    container_type_check_base<is_rev_row_vector, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_rev_row_vector_st = require_all_t<
    container_type_check_base<is_rev_row_vector, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_rev_row_vector_st = require_all_not_t<
    container_type_check_base<is_rev_row_vector, scalar_type_t, TypeCheck, Check>...>;

//STAN_ADD_REQUIRE_UNARY_INNER(rev_row_vector, is_rev_row_vector, require_eigens_types);
template <typename T>
using require_vt_rev_row_vector
    = require_t<is_rev_row_vector<value_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_vt_rev_row_vector
    = require_not_t<is_rev_row_vector<value_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_vt_rev_row_vector
    = require_all_t<is_rev_row_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_vt_rev_row_vector
    = require_any_t<is_rev_row_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_vt_rev_row_vector
    = require_all_not_t<is_rev_row_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_vt_rev_row_vector
    = require_any_not_t<is_rev_row_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename T>
using require_st_rev_row_vector
    = require_t<is_rev_row_vector<scalar_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_st_rev_row_vector
    = require_not_t<is_rev_row_vector<scalar_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_st_rev_row_vector
    = require_all_t<is_rev_row_vector<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_st_rev_row_vector
    = require_any_t<is_rev_row_vector<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_st_rev_row_vector
    = require_all_not_t<is_rev_row_vector<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_st_rev_row_vector
    = require_any_not_t<is_rev_row_vector<scalar_type_t<std::decay_t<Types>>>...>;


/** \ingroup type_trait
 * Defines a static member named value which is defined to be true
 * if the type is either a type derived from `Eigen::EigenBase` with a `Scalar`
 *  type of `var_value<double>` or a `var_value<T>` where T is derived from
 * `Eigen::EigenBase`. And the type must have a compile time constant number
 *  of columns or rows equal to 1.
 */
template <typename T, typename = void>
struct is_rev_vector : std::false_type {};

//STAN_ADD_REQUIRE_UNARY(rev_vector, is_rev_vector, require_eigens_types);
template <typename T>
using require_rev_vector_t = require_t<is_rev_vector<std::decay_t<T>>>;

template <typename T>
using require_not_rev_vector_t
    = require_not_t<is_rev_vector<std::decay_t<T>>>;

template <typename... Types>
using require_all_rev_vector_t
    = require_all_t<is_rev_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_rev_vector_t
    = require_any_t<is_rev_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_rev_vector_t
    = require_all_not_t<is_rev_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_rev_vector_t
    = require_any_not_t<is_rev_vector<std::decay_t<Types>>...>;

  
//STAN_ADD_REQUIRE_CONTAINER(rev_vector, is_rev_vector, require_eigens_types);
template <template <class...> class TypeCheck, class... Check>
using require_rev_vector_vt = require_t<
    container_type_check_base<is_rev_vector, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_rev_vector_vt = require_not_t<
    container_type_check_base<is_rev_vector, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_rev_vector_vt = require_any_t<
    container_type_check_base<is_rev_vector, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_rev_vector_vt = require_any_not_t<
    container_type_check_base<is_rev_vector, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_rev_vector_vt = require_all_t<
    container_type_check_base<is_rev_vector, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_rev_vector_vt = require_all_not_t<
    container_type_check_base<is_rev_vector, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_rev_vector_st = require_t<
    container_type_check_base<is_rev_vector, scalar_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_rev_vector_st = require_not_t<
    container_type_check_base<is_rev_vector, scalar_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_rev_vector_st = require_any_t<
    container_type_check_base<is_rev_vector, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_rev_vector_st = require_any_not_t<
    container_type_check_base<is_rev_vector, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_rev_vector_st = require_all_t<
    container_type_check_base<is_rev_vector, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_rev_vector_st = require_all_not_t<
    container_type_check_base<is_rev_vector, scalar_type_t, TypeCheck, Check>...>;

//STAN_ADD_REQUIRE_UNARY_INNER(rev_vector, is_rev_vector, require_eigens_types);
template <typename T>
using require_vt_rev_vector
    = require_t<is_rev_vector<value_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_vt_rev_vector
    = require_not_t<is_rev_vector<value_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_vt_rev_vector
    = require_all_t<is_rev_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_vt_rev_vector
    = require_any_t<is_rev_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_vt_rev_vector
    = require_all_not_t<is_rev_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_vt_rev_vector
    = require_any_not_t<is_rev_vector<value_type_t<std::decay_t<Types>>>...>;

template <typename T>
using require_st_rev_vector
    = require_t<is_rev_vector<scalar_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_st_rev_vector
    = require_not_t<is_rev_vector<scalar_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_st_rev_vector
    = require_all_t<is_rev_vector<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_st_rev_vector
    = require_any_t<is_rev_vector<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_st_rev_vector
    = require_all_not_t<is_rev_vector<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_st_rev_vector
    = require_any_not_t<is_rev_vector<scalar_type_t<std::decay_t<Types>>>...>;


}  // namespace stan
#endif
