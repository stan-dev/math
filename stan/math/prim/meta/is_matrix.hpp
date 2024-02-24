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

/*! \ingroup require_eigens_types */
/*! \defgroup matrix_types matrix  */
/*! \addtogroup matrix_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_matrix */
/*! @tparam T the type to check */
template <typename T>
using require_matrix_t = require_t<is_matrix<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_matrix */
/*! @tparam T the type to check */
template <typename T>
using require_not_matrix_t = require_not_t<is_matrix<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_matrix */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_matrix_t = require_all_t<is_matrix<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_matrix */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_matrix_t = require_any_t<is_matrix<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy @ref is_matrix */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_matrix_t
    = require_all_not_t<is_matrix<std::decay_t<Types>>...>;
/*! @} */

/*! \ingroup require_eigens_types */
/*! \defgroup matrix_types matrix  */
/*! \addtogroup matrix_types */
/*! @{ */

/*! \brief Require any of the types satisfy @ref is_matrix */
/*! and any scalar type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_matrix for and whose @ref scalar_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_any_matrix_st = require_any_t<
    container_type_check_base<is_matrix, scalar_type_t, TypeCheck, Check>...>;

/*! \brief Require all of the types does not satisfy @ref is_matrix */
/*! and all scalar types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_matrix for and whose @ref scalar_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_matrix_st = require_all_t<
    container_type_check_base<is_matrix, scalar_type_t, TypeCheck, Check>...>;

/*! \brief Require none of the types satisfy @ref is_matrix */
/*! and none of the scalar types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_matrix for and whose @ref scalar_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_not_matrix_st = require_all_not_t<
    container_type_check_base<is_matrix, scalar_type_t, TypeCheck, Check>...>;
/*! @} */

}  // namespace stan

#endif
