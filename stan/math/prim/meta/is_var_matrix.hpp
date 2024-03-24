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

/*! \ingroup require_eigens_types */
/*! \defgroup var_matrix_types var_matrix  */
/*! \addtogroup var_matrix_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_var_matrix */
/*! @tparam T the type to check */
template <typename T>
using require_var_matrix_t = require_t<is_var_matrix<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_var_matrix */
/*! @tparam T the type to check */
template <typename T>
using require_not_var_matrix_t = require_not_t<is_var_matrix<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_var_matrix */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_var_matrix_t
    = require_all_t<is_var_matrix<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_var_matrix */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_var_matrix_t
    = require_any_t<is_var_matrix<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy @ref is_var_matrix */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_var_matrix_t
    = require_all_not_t<is_var_matrix<std::decay_t<Types>>...>;
/*! @} */

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

/*! \ingroup require_eigens_types */
/*! \defgroup var_col_vector_types var_col_vector  */
/*! \addtogroup var_col_vector_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_var_col_vector */
/*! @tparam T the type to check */
template <typename T>
using require_var_col_vector_t = require_t<is_var_col_vector<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_var_col_vector */
/*! @tparam T the type to check */
template <typename T>
using require_not_var_col_vector_t
    = require_not_t<is_var_col_vector<std::decay_t<T>>>;
/*! @} */

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

/*! \ingroup require_eigens_types */
/*! \defgroup var_row_vector_types var_row_vector  */
/*! \addtogroup var_row_vector_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_var_row_vector */
/*! @tparam T the type to check */
template <typename T>
using require_var_row_vector_t = require_t<is_var_row_vector<std::decay_t<T>>>;
/*! @} */

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

/*! \ingroup require_eigens_types */
/*! \defgroup var_vector_types var_vector  */
/*! \addtogroup var_vector_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_var_vector */
/*! @tparam T the type to check */
template <typename T>
using require_var_vector_t = require_t<is_var_vector<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_var_vector */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_var_vector_t
    = require_all_t<is_var_vector<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_var_vector */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_var_vector_t
    = require_any_t<is_var_vector<std::decay_t<Types>>...>;
/*! @} */

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
