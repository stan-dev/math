#ifndef STAN_MATH_PRIM_META_IS_EIGEN_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <type_traits>

namespace stan {

/** \ingroup type_trait
 * Base implimentation to check whether a type is derived from EigenBase
 */
template <typename T, typename = void>
struct is_eigen : std::false_type {};

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

namespace internal {
template <typename T>
struct is_eigen_array_impl : std::false_type {};
template <typename T, int R, int C>
struct is_eigen_array_impl<Eigen::Array<T, R, C>> : std::true_type {};
}  // namespace internal

template <typename T>
struct is_eigen_array : internal::is_eigen_array_impl<std::decay_t<T>> {};

template <typename T>
using is_eigen_matrix_or_array
    = math::disjunction<is_eigen_matrix<T>, is_eigen_array<T>>;

namespace internal {
template <typename T>
struct is_eigen_contiguous_map_impl : std::false_type {};
template <typename T, int Opts>
struct is_eigen_contiguous_map_impl<Eigen::Map<T, Opts, Eigen::Stride<0, 0>>>
    : std::true_type {};

}  // namespace internal

template <typename T>
struct is_eigen_contiguous_map
    : internal::is_eigen_contiguous_map_impl<std::decay_t<T>> {};

}  // namespace stan
#endif
