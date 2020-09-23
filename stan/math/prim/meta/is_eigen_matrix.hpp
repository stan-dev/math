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

STAN_ADD_REQUIRE_UNARY(eigen_matrix_dynamic, is_eigen_matrix_dynamic,
                       require_eigens_types);
STAN_ADD_REQUIRE_CONTAINER(eigen_matrix_dynamic, is_eigen_matrix_dynamic,
                           require_eigens_types);

}  // namespace stan

#endif
