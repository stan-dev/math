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

/*! \ingroup require_eigens_types */
/*! \defgroup rev_matrix_types rev_matrix  */
/*! \addtogroup rev_matrix_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_rev_matrix */
/*! @tparam T the type to check */
template <typename T>
using require_rev_matrix_t = require_t<is_rev_matrix<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_rev_matrix */
/*! @tparam T the type to check */
template <typename T>
using require_not_rev_matrix_t = require_not_t<is_rev_matrix<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_rev_matrix */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_rev_matrix_t
    = require_all_t<is_rev_matrix<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_rev_matrix */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_rev_matrix_t
    = require_any_t<is_rev_matrix<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy @ref is_rev_matrix */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_rev_matrix_t
    = require_all_not_t<is_rev_matrix<std::decay_t<Types>>...>;
/*! @} */

/** \ingroup type_trait
 * Defines a static member named value which is defined to be true
 * if the type is either derived from `Eigen::EigenBase` with a `Scalar`
 *  type of `var_value<double>` or a `var_value<T>` where T is derived from
 * `Eigen::EigenBase`. And the type must have a compile time constant number
 *  of columns equal to 1.
 */
template <typename T, typename = void>
struct is_rev_col_vector : std::false_type {};

/*! \ingroup require_eigens_types */
/*! \defgroup rev_col_vector_types rev_col_vector  */
/*! \addtogroup rev_col_vector_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_rev_col_vector */
/*! @tparam T the type to check */
template <typename T>
using require_rev_col_vector_t = require_t<is_rev_col_vector<std::decay_t<T>>>;
/*! @} */

/** \ingroup type_trait
 * Defines a static member named value which is defined to be true
 * if the type is either a type derived from `Eigen::EigenBase` with a `Scalar`
 *  type of `var_value<double>` or a `var_value<T>` where T is derived from
 * `Eigen::EigenBase`. And the type must have a compile time constant number
 *  of rows equal to 1.
 */
template <typename T, typename = void>
struct is_rev_row_vector : std::false_type {};

/** \ingroup type_trait
 * Defines a static member named value which is defined to be true
 * if the type is either a type derived from `Eigen::EigenBase` with a `Scalar`
 *  type of `var_value<double>` or a `var_value<T>` where T is derived from
 * `Eigen::EigenBase`. And the type must have a compile time constant number
 *  of columns or rows equal to 1.
 */
template <typename T, typename = void>
struct is_rev_vector : std::false_type {};

/*! \ingroup require_eigens_types */
/*! \defgroup rev_vector_types rev_vector  */
/*! \addtogroup rev_vector_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_rev_vector */
/*! @tparam T the type to check */
template <typename T>
using require_rev_vector_t = require_t<is_rev_vector<std::decay_t<T>>>;
/*! @} */

}  // namespace stan
#endif
