#ifndef STAN_MATH_PRIM_SCAL_META_IS_VECTOR_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_VECTOR_HPP

#include <stan/math/prim/scal/meta/is_eigen.hpp>

namespace stan {

template <typename T>
struct is_std_vector : std::false_type {
  typedef T type;
};

namespace internal {

template <typename T, bool = is_eigen<T>::value>
struct is_eigen_col_vector_impl
    : std::integral_constant<bool, T::RowsAtCompileTime == 1> {};

template <typename T>
struct is_eigen_col_vector_impl<T, false> : std::false_type {};

template <typename T, bool = is_eigen<T>::value>
struct is_eigen_row_vector_impl
    : std::integral_constant<bool, T::ColsAtCompileTime == 1> {};

template <typename T>
struct is_eigen_row_vector_impl<T, false> : std::false_type {};
}  // namespace internal

template <typename T>
struct is_eigen_col_vector : internal::is_eigen_col_vector_impl<T> {};

template <typename T>
struct is_eigen_row_vector : internal::is_eigen_row_vector_impl<T> {};

template <typename T>
struct is_eigen_vector
    : std::integral_constant<bool, is_eigen_col_vector<T>::value
                                       || is_eigen_row_vector<T>::value> {};

template <typename T, typename = void>
struct is_vector : std::false_type {};
}  // namespace stan
#endif
