#ifndef STAN_MATH_PRIM_META_IS_EIGEN_MATRIX_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_MATRIX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <type_traits>

namespace stan {

namespace internal {

/*
 * Underlying implimenation to check if a type is derived from MatrixBase
 */
template <typename T>
struct is_eigen_matrix_impl
    : bool_constant<!(T::RowsAtCompileTime == 1 || T::ColsAtCompileTime == 1)> {
};

}  // namespace internal

/*
 * Checks whether type T is derived from Eigen::MatrixBase and has columns and
 * rows not equal to 1. If true this will have a
 * static member function named value with a type of true, else value is false.
 */
template <typename T, typename Enable = void>
struct is_eigen_matrix : std::false_type {};

template <typename T>
struct is_eigen_matrix<
    T, std::enable_if_t<std::is_base_of<
           Eigen::MatrixBase<typename std::decay_t<T>::PlainObject>,
           typename std::decay_t<T>::PlainObject>::value>>
    : bool_constant<internal::is_eigen_matrix_impl<std::decay_t<T>>::value> {};

template <typename T>
struct is_eigen_matrix<
    T, std::enable_if_t<std::is_base_of<
           Eigen::MatrixBase<typename std::decay_t<T>::MatrixType>,
           typename std::decay_t<T>::MatrixType>::value>>
    : bool_constant<internal::is_eigen_matrix_impl<std::decay_t<T>>::value> {};

}  // namespace stan

#endif
