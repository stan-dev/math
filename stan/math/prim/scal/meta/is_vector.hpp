#ifndef STAN_MATH_PRIM_SCAL_META_IS_VECTOR_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_VECTOR_HPP

#include <stan/math/prim/scal/meta/bool_constant.hpp>
#include <stan/math/prim/scal/meta/is_eigen.hpp>
#include <type_traits>

namespace stan {

/**
 * Base implimentation for checking if type is std vector
 */
template <typename T, typename = void>
struct is_std_vector : std::false_type {};

namespace internal {
/**
 * Underlying implimentation for detecting if an Eigen Matrix is a column
 * vector.
 */
template <typename T, bool = is_eigen<T>::value>
struct is_eigen_col_vector_impl
    : bool_constant<std::decay_t<T>::RowsAtCompileTime == 1> {};

/**
 * Specialization for when type is not an eigen vector.
 */
template <typename T>
struct is_eigen_col_vector_impl<T, false> : std::false_type {};

/**
 * Underlying implimentation for detecting if an Eigen Matrix is a row vector.
 */
template <typename T, bool = is_eigen<T>::value>
struct is_eigen_row_vector_impl
    : std::integral_constant<bool, std::decay_t<T>::ColsAtCompileTime == 1> {};

/**
 * Specialization for when type is not an eigen vector.
 */
template <typename T>
struct is_eigen_row_vector_impl<T, false> : std::false_type {};
}  // namespace internal

/**
 * If the input type T is an eigen matrix with 1 row at compile time this
 * has a static member with a value of true. Else this has a static
 * member with a value of false.
 */
template <typename T>
struct is_eigen_col_vector : internal::is_eigen_col_vector_impl<T> {};

/**
 * If the input type T is an eigen matrix with 1 column at compile time this
 * has a static member with a value of true. Else this has a static
 * member with a value of false.
 */
template <typename T>
struct is_eigen_row_vector : internal::is_eigen_row_vector_impl<T> {};

/**
 * If the input type T is an eigen matrix with 1 column or 1 row at compile time
 * this has a static member with a value of true. Else this has a static
 * member with a value of false.
 */
template <typename T>
struct is_eigen_vector : bool_constant<is_eigen_col_vector<T>::value
                                       || is_eigen_row_vector<T>::value> {};

/**
 * If the input type T is either an eigen matrix with 1 column or 1 row at
 * compile time or a standard vector, this has a static member with a value
 * of true. Else this has a static member with a value of false.
 */
template <typename T>
struct is_vector
    : bool_constant<is_eigen_vector<T>::value || is_std_vector<T>::value> {};

}  // namespace stan
#endif
