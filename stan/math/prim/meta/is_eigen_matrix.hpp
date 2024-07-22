#ifndef STAN_MATH_PRIM_META_IS_EIGEN_MATRIX_DYNAMIC_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_MATRIX_DYNAMIC_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/conjunction.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/is_base_pointer_convertible.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {

namespace internal {

/**
 * Underlying implimenation to check if an Eigen matrix has rows or cols not
 *  equal to 1.
 */
template <typename T, bool>
struct is_eigen_matrix_dynamic_impl : std::false_type {};

template <typename T>
struct is_eigen_matrix_dynamic_impl<T, false> : std::false_type {};

template <typename T>
struct is_eigen_matrix_dynamic_impl<T, true>
    : bool_constant<!(T::RowsAtCompileTime == 1 || T::ColsAtCompileTime == 1)> {
};

}  // namespace internal

/**
 * Checks whether type T is derived from Eigen::MatrixBase and has columns and
 * rows not equal to 1. If true this will have a
 * static member function named value with a type of true, else value is false.
 * @tparam T Type to check if it is derived from `MatrixBase` and has more than
 * 1 compile time row and column.
 * @ingroup type_trait
 */
template <typename T>
struct is_eigen_matrix_dynamic
    : bool_constant<internal::is_eigen_matrix_dynamic_impl<
          std::decay_t<T>,
          is_base_pointer_convertible<Eigen::MatrixBase, T>::value>::value> {};

/*! \ingroup require_eigens_types */
/*! \defgroup eigen_matrix_dynamic_types eigen_matrix_dynamic  */
/*! \addtogroup eigen_matrix_dynamic_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_eigen_matrix_dynamic */
/*! @tparam T the type to check */
template <typename T>
using require_eigen_matrix_dynamic_t
    = require_t<is_eigen_matrix_dynamic<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_eigen_matrix_dynamic */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_eigen_matrix_dynamic_t
    = require_all_t<is_eigen_matrix_dynamic<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_eigen_matrix_dynamic */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_eigen_matrix_dynamic_t
    = require_any_t<is_eigen_matrix_dynamic<std::decay_t<Types>>...>;
/*! @} */

/*! \ingroup require_eigens_types */
/*! \defgroup eigen_matrix_dynamic_types eigen_matrix_dynamic  */
/*! \addtogroup eigen_matrix_dynamic_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_eigen_matrix_dynamic */
/*! and value type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_eigen_matrix_dynamic for and whose
 * @ref value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_eigen_matrix_dynamic_vt
    = require_t<container_type_check_base<is_eigen_matrix_dynamic, value_type_t,
                                          TypeCheck, Check...>>;
/*! @} */

}  // namespace stan

#endif
