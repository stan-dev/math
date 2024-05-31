#ifndef STAN_MATH_PRIM_META_IS_EIGEN_MATRIX_BASE_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_MATRIX_BASE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_base_pointer_convertible.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {

/**
 * Checks whether type T is derived from Eigen::MatrixBase.
 * If true this will have a static member function named value with a type
 * of true, else value is false.
 * @tparam T Type to check if it is derived from `MatrixBase`
 * @tparam Enable used for SFINAE deduction.
 * @ingroup type_trait
 */
template <typename T>
struct is_eigen_matrix_base
    : bool_constant<is_base_pointer_convertible<Eigen::MatrixBase, T>::value> {
};

/*! \ingroup require_eigens_types */
/*! \defgroup eigen_matrix_base_types eigen_matrix_base  */
/*! \addtogroup eigen_matrix_base_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_eigen_matrix_base */
/*! @tparam T the type to check */
template <typename T>
using require_eigen_matrix_base_t
    = require_t<is_eigen_matrix_base<std::decay_t<T>>>;

/*! \brief Require type satisfies @ref is_eigen_matrix_base */
/*! and value type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_eigen_matrix_base for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_eigen_matrix_base_vt
    = require_t<container_type_check_base<is_eigen_matrix_base, value_type_t,
                                          TypeCheck, Check...>>;

/*! \brief Require all of the types satisfy is_eigen_matrix_base */
/*! and all of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_eigen_matrix_base for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_eigen_matrix_base_vt
    = require_all_t<container_type_check_base<
        is_eigen_matrix_base, value_type_t, TypeCheck, Check>...>;
/*! @} */

}  // namespace stan

#endif
