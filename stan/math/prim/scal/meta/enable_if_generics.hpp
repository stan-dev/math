#ifndef STAN_MATH_PRIM_SCAL_META_ENABLE_IF_GENERICS_HPP
#define STAN_MATH_PRIM_SCAL_META_ENABLE_IF_GENERICS_HPP

#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>
#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <stan/math/prim/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/is_vector.hpp>
#include <stan/math/prim/scal/meta/is_vector_like.hpp>
#include <stan/math/prim/scal/meta/is_eigen.hpp>
#include <stan/math/prim/scal/meta/is_fvar.hpp>
#include <stan/math/prim/scal/meta/is_var_or_arithmetic.hpp>

#include <type_traits>

namespace stan {

/**
 * Template metaprogram for deducing if a method should be remove a
 * a function from overload resolution. This method performs a negation of
 * the input boolean.
 */
template <typename Check>
using disable_if = typename std::enable_if_t<!Check::value>;

/**
 * Template metaprogram for deducing if a method should be remove a
 * a function from overload resolution.
 *
 * Returns a type void if all conditions are true and otherwise fails.
 */
template <class... Checks>
using enable_if_all = std::enable_if_t<math::conjunction<Checks...>::value>;

/**
 * Template metaprogram for deducing if a method should be remove a
 * a function from overload resolution.
 *
 * Returns a type void if any of the conditions are true and otherwise fails.
 */
template <class... Checks>
using enable_if_any = std::enable_if_t<math::disjunction<Checks...>::value>;

/**
 * Template metaprogram for deducing if a method should be remove a
 * a function from overload resolution.
 *
 * Returns a type void if all of the conditions are false.
 */
template <class... Checks>
using disable_if_all = std::enable_if_t<!math::disjunction<Checks...>::value>;

/**
 * Template metaprogram for deducing if a method should be remove a
 * a function from overload resolution.
 *
 * Returns a type void if any of the conditions are false.
 */
template <class... Checks>
using disable_if_any = std::enable_if_t<!math::conjunction<Checks...>::value>;

// Check whether the decayed types are the same
template <typename T, typename S>
struct is_same_decay
    : std::integral_constant<
          bool, std::is_same<std::decay_t<T>, std::decay_t<S>>::value> {};

template <typename T, typename S>
using enable_if_same = std::enable_if_t<is_same_decay<T, S>::value>;

template <typename T, typename S>
using disable_if_same = disable_if<is_same_decay<T, S>>;

template <typename T, typename... Types>
using enable_if_all_same = enable_if_all<is_same_decay<T, Types>...>;

template <typename T, typename... Types>
using disable_if_all_same = disable_if_all<is_same_decay<T, Types>...>;

// Checks decayed (non-const/ref'd) type is arithmetic
template <typename T>
struct is_arithmetic_decay
    : std::integral_constant<bool, std::is_arithmetic<std::decay_t<T>>::value> {
};

template <typename T>
using enable_if_arithmetic = std::enable_if_t<is_arithmetic_decay<T>::value>;

template <typename T>
using disable_if_arithmetic = disable_if<is_arithmetic_decay<T>>;

template <typename... Types>
using enable_if_all_arithmetic = enable_if_all<is_arithmetic_decay<Types>...>;

template <typename... Types>
using enable_if_any_arithmetic = enable_if_any<is_arithmetic_decay<Types>...>;

template <typename... Types>
using disable_if_all_arithmetic
    = disable_if_all<is_arithmetic_decay<Types>...>;

template <typename... Types>
using disable_if_any_arithmetic
    = disable_if_any<is_arithmetic_decay<Types>...>;

// Check if a type contains or is arithmetic
template <typename T>
struct is_contains_arithmetic
    : std::integral_constant<
          bool, std::is_arithmetic<scalar_type_t<std::decay_t<T>>>::value> {};

template <typename T>
using enable_if_contains_arithmetic
    = std::enable_if_t<is_contains_arithmetic<T>::value>;

template <typename T>
using disable_if_contains_arithmetic
    = disable_if<is_contains_arithmetic<T>>;

template <typename... Types>
using enable_if_all_contains_arithmetic
    = enable_if_all<is_contains_arithmetic<Types>...>;

template <typename... Types>
using enable_if_any_contains_arithmetic
    = enable_if_any<is_contains_arithmetic<Types>...>;

template <typename... Types>
using disable_if_all_contains_arithmetic
    = disable_if_all<is_contains_arithmetic<Types>...>;

template <typename... Types>
using disable_if_any_contains_arithmetic
    = disable_if_any<is_contains_arithmetic<Types>...>;

// Checks whether the type is floating_point
template <typename T>
struct is_fp_decay
    : std::integral_constant<bool,
                             std::is_floating_point<std::decay_t<T>>::value> {};

template <typename T>
using enable_if_floating_point = std::enable_if_t<is_fp_decay<T>::value>;

template <typename T>
using disable_if_floating_point = disable_if<is_fp_decay<T>>;

template <typename... Types>
using enable_if_all_floating_point = enable_if_all<is_fp_decay<Types>...>;

template <typename... Types>
using enable_if_any_floating_point = enable_if_any<is_fp_decay<Types>...>;

template <typename... Types>
using disable_if_all_floating_point
    = disable_if_all<is_fp_decay<Types>...>;

template <typename... Types>
using disable_if_any_floating_point
    = disable_if_any<is_fp_decay<Types>...>;

// Check if a type contains or is floating point
template <typename T>
struct is_contains_floating_point
    : std::integral_constant<
          bool, std::is_floating_point<scalar_type_t<std::decay_t<T>>>::value> {
};

template <typename T>
using enable_if_contains_floating_point
    = std::enable_if_t<is_contains_floating_point<T>::value>;

template <typename T>
using disable_if_contains_floating_point
    = disable_if<is_contains_floating_point<T>>;

template <typename... Types>
using enable_if_all_contains_floating_point
    = enable_if_all<is_contains_floating_point<Types>...>;

template <typename... Types>
using enable_if_any_contains_floating_point
    = enable_if_any<is_contains_floating_point<Types>...>;

template <typename... Types>
using disable_if_all_contains_floating_point
    = disable_if_all<is_contains_floating_point<Types>...>;

template <typename... Types>
using disable_if_any_contains_floating_point
    = disable_if_any<is_contains_floating_point<Types>...>;

// Check if type is a var

template <typename T>
using enable_if_var = std::enable_if_t<is_var<std::decay<T>>::value>;

template <typename T>
using disable_if_var = disable_if<is_var<std::decay<T>>>;

template <typename... Types>
using enable_if_all_var = enable_if_all<is_var<std::decay<Types>>...>;

template <typename... Types>
using enable_if_any_var = enable_if_any<is_var<std::decay<Types>>...>;

template <typename... Types>
using disable_if_all_var = disable_if_all<is_var<std::decay<Types>>...>;

template <typename... Types>
using disable_if_any_var = disable_if_any<is_var<std::decay<Types>>...>;

// Check if type contains or is a var
template <typename T>
struct is_contains_var
    : std::integral_constant<bool,
                             is_var<scalar_type_t<std::decay_t<T>>>::value> {};

template <typename T>
using enable_if_contains_var = std::enable_if_t<is_contains_var<T>::value>;

template <typename T>
using disable_if_contains_var = disable_if<is_contains_var<T>>;

template <typename... Types>
using enable_if_all_contains_var = enable_if_all<is_contains_var<Types>...>;

template <typename... Types>
using enable_if_any_contains_var = enable_if_any<is_contains_var<Types>...>;

template <typename... Types>
using disable_if_all_contains_var
    = disable_if_all<is_contains_var<Types>...>;

template <typename... Types>
using disable_if_any_contains_var
    = disable_if_any<is_contains_var<Types>...>;

// Check if type is a fvar
template <typename T>
using enable_if_fvar = std::enable_if_t<is_fvar<std::decay<T>>::value>;

template <typename T>
using disable_if_fvar = disable_if<is_fvar<std::decay<T>>>;

template <typename... Types>
using enable_if_all_fvar = enable_if_all<is_fvar<std::decay<Types>>...>;

template <typename... Types>
using enable_if_any_fvar = enable_if_any<is_fvar<std::decay<Types>>...>;

template <typename... Types>
using disable_if_all_fvar = disable_if_all<is_fvar<std::decay<Types>>...>;

template <typename... Types>
using disable_if_any_fvar = disable_if_any<is_fvar<std::decay<Types>>...>;

// Check if type contains or is a fvar
template <typename T>
struct is_contains_fvar
    : std::integral_constant<bool,
                             is_fvar<scalar_type_t<std::decay_t<T>>>::value> {};

template <typename T>
using enable_if_contains_fvar = std::enable_if_t<is_contains_fvar<T>::value>;

template <typename T>
using disable_if_contains_fvar = disable_if<is_contains_fvar<T>>;

template <typename... Types>
using enable_if_all_contains_fvar = enable_if_all<is_contains_fvar<Types>...>;

template <typename... Types>
using enable_if_any_contains_fvar = enable_if_any<is_contains_fvar<Types>...>;

template <typename... Types>
using disable_if_all_contains_fvar
    = disable_if_all<is_contains_fvar<Types>...>;

template <typename... Types>
using disable_if_any_contains_fvar
    = disable_if_any<is_contains_fvar<Types>...>;

template <typename T>
struct is_ad_type
    : std::integral_constant<bool, is_fvar<std::decay_t<T>>::value
                                       || is_var<std::decay_t<T>>::value> {};

template <typename T>
using enable_if_ad_type = std::enable_if_t<is_ad_type<std::decay<T>>::value>;

template <typename T>
using disable_if_ad_type = disable_if<is_ad_type<std::decay<T>>>;

template <typename... Types>
using enable_if_all_ad_type = enable_if_all<is_ad_type<std::decay<Types>>...>;

template <typename... Types>
using enable_if_any_ad_type = enable_if_any<is_ad_type<std::decay<Types>>...>;

template <typename... Types>
using disable_if_all_ad_type
    = disable_if_all<is_ad_type<std::decay<Types>>...>;

template <typename... Types>
using disable_if_any_ad_type
    = disable_if_any<is_ad_type<std::decay<Types>>...>;

// Check if type contains or is a ad_type
template <typename T>
struct is_contains_ad_type
    : std::integral_constant<
          bool, is_ad_type<scalar_type_t<std::decay_t<T>>>::value> {};

template <typename T>
using enable_if_contains_ad_type
    = std::enable_if_t<is_contains_ad_type<T>::value>;

template <typename T>
using disable_if_contains_ad_type = disable_if<is_contains_ad_type<T>>;

template <typename... Types>
using enable_if_all_contains_ad_type
    = enable_if_all<is_contains_ad_type<Types>...>;

template <typename... Types>
using enable_if_any_contains_ad_type
    = enable_if_any<is_contains_ad_type<Types>...>;

template <typename... Types>
using disable_if_all_contains_ad_type
    = disable_if_all<is_contains_ad_type<Types>...>;

template <typename... Types>
using disable_if_any_contains_ad_type
    = disable_if_any<is_contains_ad_type<Types>...>;

// Enables if type is var or arithmetic
template <typename T>
using enable_if_var_or_arithmetic
    = std::enable_if_t<is_var_or_arithmetic<T>::value>;

template <typename T>
using disable_if_var_or_arithmetic = disable_if<is_var_or_arithmetic<T>>;

template <typename... Types>
using enable_if_all_var_or_arithmetic
    = enable_if_all<is_var_or_arithmetic<Types>...>;

template <typename... Types>
using enable_if_any_var_or_arithmetic
    = enable_if_any<is_var_or_arithmetic<Types>...>;

template <typename... Types>
using disable_if_all_var_or_arithmetic
    = disable_if_all<is_var_or_arithmetic<Types>...>;

template <typename... Types>
using disable_if_any_var_or_arithmetic
    = disable_if_any<is_var_or_arithmetic<Types>...>;

// Enables if type is var or arithmetic
template <typename T>
struct is_contains_var_or_arithmetic
    : std::integral_constant<bool,
                             is_var_or_arithmetic<scalar_type_t<T>>::value> {};

template <typename T>
using enable_if_contains_var_or_arithmetic
    = std::enable_if_t<is_contains_var_or_arithmetic<T>::value>;

template <typename T>
using disable_if_contains_var_or_arithmetic
    = disable_if<is_contains_var_or_arithmetic<T>>;

template <typename... Types>
using enable_if_all_contains_var_or_arithmetic
    = enable_if_all<is_contains_var_or_arithmetic<Types>...>;

template <typename... Types>
using enable_if_any_contains_var_or_arithmetic
    = enable_if_any<is_contains_var_or_arithmetic<Types>...>;

template <typename... Types>
using disable_if_all_contains_var_or_arithmetic
    = disable_if_all<is_contains_var_or_arithmetic<Types>...>;

template <typename... Types>
using disable_if_any_contains_var_or_arithmetic
    = disable_if_any<is_contains_var_or_arithmetic<Types>...>;

// Checks whether type is arithmetic, var, or fvar
template <typename T>
struct is_stan_scalar
    : std::integral_constant<bool, std::is_arithmetic<std::decay_t<T>>::value
                                       || is_var<std::decay_t<T>>::value
                                       || is_fvar<std::decay_t<T>>::value> {};

template <typename T>
using enable_if_stan_scalar = std::enable_if_t<is_stan_scalar<T>::value>;

template <typename T>
using disable_if_stan_scalar = disable_if<is_stan_scalar<T>>;

template <typename... Types>
using enable_if_all_stan_scalar = enable_if_all<is_stan_scalar<Types>...>;

template <typename... Types>
using enable_if_any_stan_scalar = enable_if_any<is_stan_scalar<Types>...>;

template <typename... Types>
using disable_if_all_stan_scalar
    = disable_if_all<is_stan_scalar<Types>...>;

template <typename... Types>
using disable_if_any_stan_stan_scalar
    = disable_if_any<is_stan_scalar<Types>...>;

// Check whether a type contains (or is) arithmetic, var, or fvar
template <typename T>
struct is_contains_stan_scalar
    : std::integral_constant<bool, is_stan_scalar<scalar_type_t<T>>::value> {};

template <typename T>
using enable_if_contains_stan_scalar
    = std::enable_if_t<is_contains_stan_scalar<T>::value>;

template <typename T>
using disable_if_contains_stan_scalar
    = disable_if<is_contains_stan_scalar<T>>;

template <typename... Types>
using enable_if_all_contains_stan_scalar
    = enable_if_all<is_contains_stan_scalar<Types>...>;

template <typename... Types>
using enable_if_any_contains_stan_scalar
    = enable_if_any<is_contains_stan_scalar<Types>...>;

template <typename... Types>
using disable_if_all_contains_stan_scalar
    = enable_if_all<is_contains_stan_scalar<Types>...>;

template <typename... Types>
using disable_if_any_stan_contains_stan_scalar
    = enable_if_any<is_contains_stan_scalar<Types>...>;

// Checks whether type is a scalar as defined by the standard
template <typename T>
struct is_scalar_decay
    : std::integral_constant<bool, std::is_scalar<std::decay_t<T>>::value> {};

template <typename T>
using enable_if_scalar = std::enable_if_t<is_scalar_decay<T>::value>;

template <typename T>
using disable_if_scalar = disable_if<is_scalar_decay<T>>;

template <typename... Types>
using enable_if_all_scalar = enable_if_all<is_scalar_decay<Types>...>;

template <typename... Types>
using enable_if_any_scalar = enable_if_any<is_scalar_decay<Types>...>;

template <typename... Types>
using disable_if_all_scalar = disable_if_all<is_scalar_decay<Types>...>;

template <typename... Types>
using disable_if_any_scalar = disable_if_any<is_scalar_decay<Types>...>;


// Checks whether decayed type is a vector
template <typename T>
struct is_vector_decay
    : std::integral_constant<bool, is_vector<std::decay_t<T>>::value> {};

template <typename T>
using enable_if_vector = std::enable_if_t<is_vector_decay<T>::value>;

template <typename T>
using disable_if_vector = disable_if<is_vector_decay<T>>;

template <typename... Types>
using enable_if_all_vector = enable_if_all<is_vector_decay<Types>...>;

template <typename... Types>
using enable_if_any_vector = enable_if_any<is_vector_decay<Types>...>;

template <typename... Types>
using disable_if_all_vector = disable_if_all<is_vector_decay<Types>...>;

template <typename... Types>
using disable_if_any_vector = disable_if_any<is_vector_decay<Types>...>;

// Checks whether decayed type is a standard vector
template <typename T>
struct is_std_vector_decay
    : std::integral_constant<bool, is_std_vector<std::decay_t<T>>::value> {};

template <typename T>
using enable_if_std_vector = std::enable_if_t<is_std_vector_decay<T>::value>;

template <typename T>
using disable_if_std_vector = disable_if<is_std_vector_decay<T>>;

template <typename... Types>
using enable_if_all_std_vector = enable_if_all<is_std_vector_decay<Types>...>;

template <typename... Types>
using enable_if_any_std_vector = enable_if_any<is_std_vector_decay<Types>...>;

template <typename... Types>
using disable_if_all_std_vector
    = disable_if_all<is_std_vector_decay<Types>...>;

template <typename... Types>
using disable_if_any_std_vector
    = disable_if_any<is_std_vector_decay<Types>...>;

// Enable if function comes from Eigen
template <typename T>
using enable_if_eigen = std::enable_if_t<is_eigen<T>::value>;

template <typename T>
using disable_if_eigen = disable_if<is_eigen<T>>;

template <typename... Types>
using enable_if_all_eigen = enable_if_all<is_eigen<Types>...>;

template <typename... Types>
using enable_if_any_eigen = enable_if_any<is_eigen<Types>...>;

template <typename... Types>
using disable_if_all_eigen = disable_if_all<is_eigen<Types>...>;

template <typename... Types>
using disable_if_any_eigen = disable_if_any<is_eigen<Types>...>;

// Checks if type is Eigen or arithmetic, var, or fvar
template <typename T>
struct is_eigen_or_stan_scalar
    : std::integral_constant<bool, is_eigen<T>::value
                                       || is_stan_scalar<T>::value> {};

template <typename T>
using enable_if_eigen_or_stan_scalar
    = std::enable_if_t<is_eigen_or_stan_scalar<T>::value>;

template <typename T>
using disable_if_eigen_or_stan_scalar
    = disable_if<is_eigen_or_stan_scalar<T>>;

template <typename... Types>
using enable_if_all_eigen_or_stan_scalar
    = enable_if_all<is_eigen_or_stan_scalar<Types>...>;

template <typename... Types>
using enable_if_any_eigen_or_stan_scalar
    = enable_if_any<is_eigen_or_stan_scalar<Types>...>;

template <typename... Types>
using disable_if_all_eigen_or_stan_scalar
    = disable_if_all<is_eigen_or_stan_scalar<Types>...>;

template <typename... Types>
using disable_if_any_eigen_or_stan_scalar
    = disable_if_any<is_eigen_or_stan_scalar<Types>...>;

// Checks if type is Eigen scalar type and arithmetic
template <typename T>
struct is_eigen_arithmetic
    : std::integral_constant<bool, is_eigen<T>::value
                                       && is_contains_arithmetic<T>::value> {};

template <typename T>
using enable_if_eigen_arithmetic
    = std::enable_if_t<is_eigen_arithmetic<T>::value>;

template <typename T>
using disable_if_eigen_arithmetic = disable_if<is_eigen_arithmetic<T>>;

template <typename... Types>
using enable_if_all_eigen_arithmetic
    = enable_if_all<is_eigen_arithmetic<Types>...>;

template <typename... Types>
using enable_if_any_eigen_arithmetic
    = enable_if_any<is_eigen_arithmetic<Types>...>;

template <typename... Types>
using disable_if_all_eigen_arithmetic
    = disable_if_all<is_eigen_arithmetic<Types>...>;

template <typename... Types>
using disable_if_any_eigen_arithmetic
    = disable_if_any<is_eigen_arithmetic<Types>...>;

// Enable if for Eigen col vectors
template <typename T>
using enable_if_eigen_col_vector
    = std::enable_if_t<is_eigen_col_vector<T>::value>;

template <typename T>
using disable_if_eigen_col_vector = disable_if<is_eigen_col_vector<T>>;

template <typename... Types>
using enable_if_all_eigen_col_vector
    = enable_if_all<is_eigen_col_vector<Types>...>;

template <typename... Types>
using enable_if_any_eigen_col_vector
    = enable_if_any<is_eigen_col_vector<Types>...>;

template <typename... Types>
using disable_if_all_eigen_col_vector
    = disable_if_all<is_eigen_col_vector<Types>...>;

template <typename... Types>
using disable_if_any_eigen_col_vector
    = disable_if_any<is_eigen_col_vector<Types>...>;

// Enable if for Eigen row vectors
template <typename T>
using enable_if_eigen_row_vector
    = std::enable_if_t<is_eigen_row_vector<T>::value>;

template <typename T>
using disable_if_eigen_row_vector = disable_if<is_eigen_row_vector<T>>;

template <typename... Types>
using enable_if_all_eigen_row_vector
    = enable_if_all<is_eigen_row_vector<Types>...>;

template <typename... Types>
using enable_if_any_eigen_row_vector
    = enable_if_any<is_eigen_row_vector<Types>...>;

template <typename... Types>
using disable_if_all_eigen_row_vector
    = disable_if_all<is_eigen_row_vector<Types>...>;

template <typename... Types>
using disable_if_any_eigen_row_vector
    = disable_if_any<is_eigen_row_vector<Types>...>;

// Enable if eigen row or column vector
template <typename T>
using enable_if_eigen_vector = std::enable_if_t<is_eigen_vector<T>::value>;

template <typename T>
using disable_if_eigen_vector = disable_if<is_eigen_vector<T>>;

template <typename... Types>
using enable_if_all_eigen_vector = enable_if_all<is_eigen_vector<Types>...>;

template <typename... Types>
using enable_if_any_eigen_vector = enable_if_any<is_eigen_vector<Types>...>;

template <typename... Types>
using disable_if_all_eigen_vector
    = disable_if_all<is_eigen_vector<Types>...>;

template <typename... Types>
using disable_if_any_eigen_vector
    = disable_if_any<is_eigen_vector<Types>...>;

// Check whether Eigen types satisfy a dot product
template <typename T1, typename T2>
struct is_dot_product
    : std::integral_constant<bool, is_eigen_row_vector<T1>::value
                                       && is_eigen_col_vector<T2>::value> {};

template <typename T1, typename T2>
using enable_if_dot_product = std::enable_if_t<is_dot_product<T1, T2>::value>;

template <typename T1, typename T2>
using disable_if_dot_product = disable_if<is_dot_product<T1, T2>>;

namespace internal {
// primary template for checking if eigen matrix rows match
template <class T1, class T2, bool = is_eigen<T1>::value,
          bool = is_eigen<T2>::value>
struct is_eigen_rows_match_impl
    : std::integral_constant<bool,
                             T1::RowsAtCompileTime == T2::RowsAtCompileTime> {};

// if not eigen
template <class T1, class T2>
struct is_eigen_rows_match_impl<T1, T2, false, false> : std::false_type {};

template <class T1, class T2>
struct is_eigen_rows_match_impl<T1, T2, true, false> : std::false_type {};

template <class T1, class T2>
struct is_eigen_rows_match_impl<T1, T2, false, true> : std::false_type {};

}  // namespace internal

// Check whether eigen matrices rows match
template <typename T1, typename T2>
struct is_eigen_rows_match : internal::is_eigen_rows_match_impl<T1, T2> {};

// Enables for matching rows and columns
template <typename T1, typename T2>
using enable_if_eigen_rows_match
    = std::enable_if_t<is_eigen_rows_match<T1, T2>::value>;

template <typename T1, typename T2>
using disable_if_eigen_rows_match
    = disable_if<is_eigen_rows_match<T1, T2>>;

namespace internal {
// primary template for checking if eigen matrix cols match
template <class T1, class T2, bool = is_eigen<T1>::value,
          bool = is_eigen<T2>::value>
struct is_eigen_cols_match_impl
    : std::integral_constant<bool,
                             T1::ColsAtCompileTime == T2::ColsAtCompileTime> {};

// if not eigen
template <class T1, class T2>
struct is_eigen_cols_match_impl<T1, T2, false, false> : std::false_type {};

template <class T1, class T2>
struct is_eigen_cols_match_impl<T1, T2, true, false> : std::false_type {};

template <class T1, class T2>
struct is_eigen_cols_match_impl<T1, T2, false, true> : std::false_type {};

}  // namespace internal

// Check whether eigen matrices cols match
template <typename T1, typename T2>
struct is_eigen_cols_match : internal::is_eigen_cols_match_impl<T1, T2> {};

template <typename T1, typename T2>
using enable_if_eigen_cols_match
    = std::enable_if_t<is_eigen_cols_match<T1, T2>::value>;

template <typename T1, typename T2>
using disable_if_eigen_cols_match
    = disable_if<is_eigen_cols_match<T1, T2>>;

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
using disable_if_eigen_dims_match
    = disable_if<is_eigen_dims_match<T1, T2>>;

}  // namespace stan
#endif
