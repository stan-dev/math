#ifndef STAN_MATH_PRIM_MAT_META_IS_VECTOR_HPP
#define STAN_MATH_PRIM_MAT_META_IS_VECTOR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/meta/is_eigen.hpp>
#include <stan/math/prim/scal/meta/is_vector.hpp>

namespace stan {

namespace internal {

// primary template for checking if eigen matrix is column vector
template <class T, typename = void>
struct is_eigen_col_vector_impl : std::false_type {};

template <class T, typename = void>
struct is_eigen_row_vector_impl : std::false_type {};

// Template that checks Eigen matrix for one column
template <class T>
struct is_eigen_col_vector_impl<
    T, std::enable_if_t<is_eigen<std::decay_t<T>>::value>>
    : std::integral_constant<bool, T::ColsAtCompileTime == 1> {};

// Template that checks Eigen matrix for one row
template <class T>
struct is_eigen_row_vector_impl<
    T, std::enable_if_t<is_eigen<std::decay_t<T>>::value>>
    : std::integral_constant<bool, T::RowsAtCompileTime == 1> {};

}  // namespace internal

// Check whether type is an eigen col vector
template <typename T>
struct is_eigen_col_vector<
    T, std::enable_if_t<internal::is_eigen_col_vector_impl<T>::value>>
    : std::true_type {};

// Check whether type is an eigen row vector
template <typename T>
struct is_eigen_row_vector<
    T, std::enable_if_t<internal::is_eigen_row_vector_impl<T>::value>>
    : std::true_type {};

// Checks whether decayed type is an eigen vector
template <typename T>
struct is_eigen_vector<T, std::enable_if_t<is_eigen_col_vector<T>::value
                                           || is_eigen_row_vector<T>::value>>
    : std::true_type {};

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
