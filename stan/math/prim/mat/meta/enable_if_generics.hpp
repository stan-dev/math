#ifndef STAN_MATH_PRIM_MAT_META_ENABLE_IF_GENERICS_HPP
#define STAN_MATH_PRIM_MAT_META_ENABLE_IF_GENERICS_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/meta/is_eigen.hpp>
#include <stan/math/prim/mat/meta/is_vector.hpp>
#include <stan/math/prim/arr/meta/is_vector.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>
#include <stan/math/prim/scal/meta/enable_if_generics.hpp>
#include <type_traits>

namespace stan {

// Checks whether decayed type is a vector
template <typename T>
struct is_vector_decay
    : std::integral_constant<bool, is_vector<std::decay_t<T>>::value> {};

template <typename T>
using enable_if_vector = std::enable_if_t<is_vector_decay<T>::value>;

template <typename T>
using enable_if_not_vector = enable_if_not<is_vector_decay<T>>;

template <typename... Types>
using enable_if_all_vector = enable_if_all<is_vector_decay<Types>...>;

template <typename... Types>
using enable_if_any_vector = enable_if_any<is_vector_decay<Types>...>;

template <typename... Types>
using enable_if_all_not_vector = enable_if_all_not<is_vector_decay<Types>...>;

template <typename... Types>
using enable_if_any_not_vector = enable_if_any_not<is_vector_decay<Types>...>;

// Checks whether decayed type is a standard vector
template <typename T>
struct is_std_vector_decay
    : std::integral_constant<bool, is_std_vector<std::decay_t<T>>::value> {};

template <typename T>
using enable_if_std_vector = std::enable_if_t<is_std_vector_decay<T>::value>;

template <typename T>
using enable_if_not_std_vector = enable_if_not<is_std_vector_decay<T>>;

template <typename... Types>
using enable_if_all_std_vector = enable_if_all<is_std_vector_decay<Types>...>;

template <typename... Types>
using enable_if_any_std_vector = enable_if_any<is_std_vector_decay<Types>...>;

template <typename... Types>
using enable_if_all_not_std_vector
    = enable_if_all_not<is_std_vector_decay<Types>...>;

template <typename... Types>
using enable_if_any_not_std_vector
    = enable_if_any_not<is_std_vector_decay<Types>...>;

// Enable if function comes from Eigen
template <typename T>
using enable_if_eigen = std::enable_if_t<is_eigen_decay<T>::value>;

template <typename T>
using enable_if_not_eigen = enable_if_not<is_eigen_decay<T>>;

template <typename... Types>
using enable_if_all_eigen = enable_if_all<is_eigen_decay<Types>...>;

template <typename... Types>
using enable_if_any_eigen = enable_if_any<is_eigen_decay<Types>...>;

template <typename... Types>
using enable_if_all_not_eigen = enable_if_all_not<is_eigen_decay<Types>...>;

template <typename... Types>
using enable_if_any_not_eigen = enable_if_any_not<is_eigen_decay<Types>...>;

// Checks if type is Eigen or arithmetic, var, or fvar
template <typename T>
struct is_eigen_or_stan_scalar
    : std::integral_constant<bool, is_eigen_decay<T>::value
                                       || is_stan_scalar<T>::value> {};

template <typename T>
using enable_if_eigen_or_stan_scalar
    = std::enable_if_t<is_eigen_or_stan_scalar<T>::value>;

template <typename T>
using enable_if_not_eigen_or_stan_scalar
    = enable_if_not<is_eigen_or_stan_scalar<T>>;

template <typename... Types>
using enable_if_all_eigen_or_stan_scalar
    = enable_if_all<is_eigen_or_stan_scalar<Types>...>;

template <typename... Types>
using enable_if_any_eigen_or_stan_scalar
    = enable_if_any<is_eigen_or_stan_scalar<Types>...>;

template <typename... Types>
using enable_if_all_not_eigen_or_stan_scalar
    = enable_if_all_not<is_eigen_or_stan_scalar<Types>...>;

template <typename... Types>
using enable_if_any_not_eigen_or_stan_scalar
    = enable_if_any_not<is_eigen_or_stan_scalar<Types>...>;

// Checks if type is Eigen scalar type is arithmetic
template <typename T>
struct is_eigen_arithmetic
    : std::integral_constant<bool, is_eigen_decay<T>::value
                                       && is_contains_arithmetic<T>::value> {};

template <typename T>
using enable_if_eigen_arithmetic
    = std::enable_if_t<is_eigen_arithmetic<T>::value>;

template <typename T>
using enable_if_not_eigen_arithmetic = enable_if_not<is_eigen_arithmetic<T>>;

template <typename... Types>
using enable_if_all_eigen_arithmetic
    = enable_if_all<is_eigen_arithmetic<Types>...>;

template <typename... Types>
using enable_if_any_eigen_arithmetic
    = enable_if_any<is_eigen_arithmetic<Types>...>;

template <typename... Types>
using enable_if_all_not_eigen_arithmetic
    = enable_if_all_not<is_eigen_arithmetic<Types>...>;

template <typename... Types>
using enable_if_any_not_eigen_arithmetic
    = enable_if_any_not<is_eigen_arithmetic<Types>...>;

// Enable if for Eigen col vectors
template <typename T>
using enable_if_eigen_col_vector
    = std::enable_if_t<is_eigen_col_vector<T>::value>;

template <typename T>
using enable_if_not_eigen_col_vector = enable_if_not<is_eigen_col_vector<T>>;

template <typename... Types>
using enable_if_all_eigen_col_vector
    = enable_if_all<is_eigen_col_vector<Types>...>;

template <typename... Types>
using enable_if_any_eigen_col_vector
    = enable_if_any<is_eigen_col_vector<Types>...>;

template <typename... Types>
using enable_if_all_not_eigen_col_vector
    = enable_if_all_not<is_eigen_col_vector<Types>...>;

template <typename... Types>
using enable_if_any_not_eigen_col_vector
    = enable_if_any_not<is_eigen_col_vector<Types>...>;

// Enable if for Eigen row vectors
template <typename T>
using enable_if_eigen_row_vector
    = std::enable_if_t<is_eigen_row_vector<T>::value>;

template <typename T>
using enable_if_not_eigen_row_vector = enable_if_not<is_eigen_row_vector<T>>;

template <typename... Types>
using enable_if_all_eigen_row_vector
    = enable_if_all<is_eigen_row_vector<Types>...>;

template <typename... Types>
using enable_if_any_eigen_row_vector
    = enable_if_any<is_eigen_row_vector<Types>...>;

template <typename... Types>
using enable_if_all_not_eigen_row_vector
    = enable_if_all_not<is_eigen_row_vector<Types>...>;

template <typename... Types>
using enable_if_any_not_eigen_row_vector
    = enable_if_any_not<is_eigen_row_vector<Types>...>;

// Enable if eigen row or column vector
template <typename T>
using enable_if_eigen_vector = std::enable_if_t<is_eigen_vector<T>::value>;

template <typename T>
using enable_if_not_eigen_vector = enable_if_not<is_eigen_vector<T>>;

template <typename... Types>
using enable_if_all_eigen_vector = enable_if_all<is_eigen_vector<Types>...>;

template <typename... Types>
using enable_if_any_eigen_vector = enable_if_any<is_eigen_vector<Types>...>;

template <typename... Types>
using enable_if_all_not_eigen_vector
    = enable_if_all_not<is_eigen_vector<Types>...>;

template <typename... Types>
using enable_if_any_not_eigen_vector
    = enable_if_any_not<is_eigen_vector<Types>...>;

// Check whether Eigen types satisfy a dot product
template <typename T1, typename T2>
struct is_dot_product
    : std::integral_constant<bool, is_eigen_row_vector<T1>::value
                                       && is_eigen_col_vector<T2>::value> {};

template <typename T1, typename T2>
using enable_if_dot_product = std::enable_if_t<is_dot_product<T1, T2>::value>;

template <typename T1, typename T2>
using enable_if_not_dot_product = enable_if_not<is_dot_product<T1, T2>>;

namespace internal {
// primary template for checking if eigen matrix rows match
template <class T1, class T2,
          bool = is_eigen_decay<T1>::value&& is_eigen_decay<T2>::value>
struct is_eigen_rows_match_impl
    : std::integral_constant<bool,
                             T1::RowsAtCompileTime == T2::RowsAtCompileTime> {};

// if not eigen
template <class T1, class T2>
struct is_eigen_rows_match_impl<T1, T2, false> : std::false_type {};

}  // namespace internal

// Check whether eigen matrices rows match
template <typename T1, typename T2>
struct is_eigen_rows_match : internal::is_eigen_rows_match_impl<T1, T2> {};

// Enables for matching rows and columns
template <typename T1, typename T2>
using enable_if_eigen_rows_match
    = std::enable_if_t<is_eigen_rows_match<T1, T2>::value>;

template <typename T1, typename T2>
using enable_if_not_eigen_rows_match
    = enable_if_not<is_eigen_rows_match<T1, T2>>;

namespace internal {
// primary template for checking if eigen matrix cols match
template <class T1, class T2,
          bool = is_eigen_decay<T1>::value&& is_eigen_decay<T2>::value>
struct is_eigen_cols_match_impl
    : std::integral_constant<bool,
                             T1::ColsAtCompileTime == T2::ColsAtCompileTime> {};

// if not eigen
template <class T1, class T2>
struct is_eigen_cols_match_impl<T1, T2, false> : std::false_type {};

}  // namespace internal

// Check whether eigen matrices cols match
template <typename T1, typename T2>
struct is_eigen_cols_match : internal::is_eigen_cols_match_impl<T1, T2> {};

template <typename T1, typename T2>
using enable_if_eigen_cols_match
    = std::enable_if_t<is_eigen_cols_match<T1, T2>::value>;

template <typename T1, typename T2>
using enable_if_not_eigen_cols_match
    = enable_if_not<is_eigen_cols_match<T1, T2>>;

// primary template for checking if eigen matrix cols match
template <class T1, class T2>
struct is_eigen_dims_match
    : std::integral_constant<bool, is_eigen_cols_match<T1, T2>::value
                                       && is_eigen_rows_match<T1, T2>::value> {
};

template <typename T1, typename T2>
using enable_if_eigen_dims_match
    = std::enable_if_t<is_eigen_dims_match<T1, T2>::value>;

template <typename T1, typename T2>
using enable_if_not_eigen_dims_match
    = enable_if_not<is_eigen_dims_match<T1, T2>>;

}  // namespace stan
#endif
