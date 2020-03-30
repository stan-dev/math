#ifndef STAN_MATH_PRIM_META_IS_EIGEN_ROW_VECTOR_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_ROW_VECTOR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_eigen_dense.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {
/** \addtogroup type_trait
 *  @{
 */

namespace internal {
/**
 * Underlying implimentation for detecting if an Eigen Matrix is a row vector.
 */
template <typename T, bool = stan::is_eigen_dense<T>::value>
struct is_eigen_row_vector_impl
    : std::integral_constant<bool, std::decay_t<T>::RowsAtCompileTime == 1> {};

/**
 * Specialization for when type is not an eigen vector.
 */
template <typename T>
struct is_eigen_row_vector_impl<T, false> : std::false_type {};

}  // namespace internal

/** \ingroup type_trait
 * If the input type T is an eigen matrix with 1 column at compile time this
 * has a static member with a value of true. Else this has a static
 * member with a value of false.
 */
template <typename T>
struct is_eigen_row_vector : internal::is_eigen_row_vector_impl<T> {};
/** @}*/

STAN_ADD_REQUIRE_UNARY(eigen_row_vector, is_eigen_row_vector,
                       require_eigens_types);
STAN_ADD_REQUIRE_CONTAINER(eigen_row_vector, is_eigen_row_vector,
                           require_eigens_types);

}  // namespace stan

#endif
