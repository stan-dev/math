#ifndef STAN_MATH_PRIM_MAT_META_IS_VECTOR_HPP
#define STAN_MATH_PRIM_MAT_META_IS_VECTOR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/meta/is_eigen.hpp>
#include <stan/math/prim/scal/meta/is_vector.hpp>

namespace stan {

namespace internal {
// primary template for checking if eigen matrix is column vector
template <class T, bool = is_eigen_decay<T>::value>
struct is_eigen_col_vector_impl
    : std::integral_constant<bool, T::ColsAtCompileTime == 1> {};

// not column vector
template <class T>
struct is_eigen_col_vector_impl<T, false> : std::false_type {};
}  // namespace internal

// Check whether type is an eigen col vector
template <typename T>
struct is_eigen_col_vector : internal::is_eigen_col_vector_impl<T> {};

namespace internal {
// primary template for checking if eigen matrix is row vector
template <class T, bool = is_eigen_decay<T>::value>
struct is_eigen_row_vector_impl
    : std::integral_constant<bool, T::RowsAtCompileTime == 1> {};

// not row vector
template <class T>
struct is_eigen_row_vector_impl<T, false> : std::false_type {};
}  // namespace internal

// Check whether type is an eigen row vector
template <typename T>
struct is_eigen_row_vector : internal::is_eigen_row_vector_impl<T> {};

// Checks whether decayed type is an eigen vector
template <typename T>
struct is_eigen_vector
    : std::integral_constant<bool, is_eigen_col_vector<T>::value
                                       || is_eigen_row_vector<T>::value> {};

/**
 * Specialization of is_vector for Eigen vectors
 */
template <typename T>
struct is_vector<T, std::enable_if_t<is_eigen_vector<T>::value>>
    : std::true_type {
  typedef std::decay_t<T> type;
};

}  // namespace stan
#endif
