#ifndef STAN_MATH_PRIM_MAT_META_IS_EIGEN_HPP
#define STAN_MATH_PRIM_MAT_META_IS_EIGEN_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/meta/is_eigen.hpp>
#include <type_traits>

namespace stan {

namespace internal {
/*
 * Underlying implimenation to check if a type is derived from EigenBase
 */
template <typename T>
struct is_eigen_base
    : std::integral_constant<bool,
                             std::is_base_of<Eigen::EigenBase<T>, T>::value> {};

template <typename T>
struct is_eigen_base<Eigen::DenseBase<T>> : std::true_type {};

template <typename T>
struct is_eigen_base<Eigen::MatrixBase<T>> : std::true_type {};

template <typename T>
struct is_eigen_base<Eigen::EigenBase<T>> : std::true_type {};

}  // namespace internal

/*
 * Checks whether type T is derived from EigenBase. If true this will have a
 * static member function named value with a type of true, else value is false.
 */
template <typename T>
struct is_eigen<
    T, std::enable_if_t<internal::is_eigen_base<std::decay_t<T>>::value>>
    : std::true_type {};

namespace internal {
template <typename T>
struct is_eigen_matrix_impl : std::false_type {};
template <typename T, int R, int C>
struct is_eigen_matrix_impl<Eigen::Matrix<T, R, C>> : std::true_type {};
template <typename T>
struct is_eigen_matrix_impl<Eigen::SparseMatrix<T>> : std::true_type {};

}  // namespace internal

template <typename T>
struct is_eigen_matrix : internal::is_eigen_matrix_impl<std::decay_t<T>> {};

}  // namespace stan
#endif
