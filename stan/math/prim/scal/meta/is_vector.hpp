#ifndef STAN_MATH_PRIM_SCAL_META_IS_VECTOR_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_VECTOR_HPP

#include <stan/math/prim/scal/meta/is_eigen.hpp>
#include <type_traits>

namespace stan {

/**
 * Metaprogram structure to determine if type is a vector
 *
 * <p>This base class should be specialized for structured types.</p>
 *
 * @tparam T Type of object.
 */
template <typename T, typename = void>
struct is_vector : std::false_type {
  typedef std::decay_t<T> type;
};

/**
 * Metaprogram structure to determine if type is a standard vector
 *
 * <p>This base class should be specialized for structured types.</p>
 *
 * @tparam T Type of object.
 */
template <typename T, typename = void>
struct is_std_vector : std::false_type {};

namespace internal {
// primary template for checking if eigen matrix is column vector
template <class T, bool = is_eigen<std::decay_t<T>>::value>
struct is_eigen_col_vector_impl
    : std::integral_constant<bool, T::ColsAtCompileTime == 1> {};

// not eigen type
template <class T>
struct is_eigen_col_vector_impl<T, false> : std::false_type {};
}  // namespace internal

// Check whether type is an eigen col vector
template <typename T>
struct is_eigen_col_vector : internal::is_eigen_col_vector_impl<T> {};

namespace internal {
// primary template for checking if eigen matrix is row vector
template <class T, bool = is_eigen<std::decay_t<T>>::value>
struct is_eigen_row_vector_impl
    : std::integral_constant<bool, T::RowsAtCompileTime == 1> {};

// not Eigen type
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

}  // namespace stan
#endif
