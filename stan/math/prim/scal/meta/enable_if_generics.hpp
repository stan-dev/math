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
using disable_if_all = std::enable_if_t<!math::conjunction<Checks...>::value>;

/**
 * Template metaprogram for deducing if a method should be remove a
 * a function from overload resolution.
 *
 * Returns a type void if any of the conditions are false.
 */
template <class... Checks>
using disable_if_any = std::enable_if_t<!math::disjunction<Checks...>::value>;

// Check whether the decayed types are the same
template <typename T, typename S>
struct is_same_decay
    : std::integral_constant<
          bool, std::is_same<std::decay_t<T>, std::decay_t<S>>::value> {};

template <typename T, typename S>
using same_type = std::enable_if_t<is_same_decay<T, S>::value>;

template <typename T, typename S>
using not_same_type = disable_if<is_same_decay<T, S>>;

template <typename T, typename... Types>
using all_same_type = enable_if_all<is_same_decay<T, Types>...>;

template <typename T, typename... Types>
using not_all_same_type = disable_if_all<is_same_decay<T, Types>...>;

// Checks decayed (non-const/ref'd) type is arithmetic
template <typename T>
struct is_arithmetic_decay
    : std::integral_constant<bool, std::is_arithmetic<std::decay_t<T>>::value> {
};

template <typename T>
using arithmetic_type = std::enable_if_t<is_arithmetic_decay<T>::value>;

template <typename T>
using not_arithmetic_type = disable_if<is_arithmetic_decay<T>>;

template <typename... Types>
using all_arithmetic_type = enable_if_all<is_arithmetic_decay<Types>...>;

template <typename... Types>
using any_arithmetic_type = enable_if_any<is_arithmetic_decay<Types>...>;

template <typename... Types>
using not_all_arithmetic_type = disable_if_all<is_arithmetic_decay<Types>...>;

template <typename... Types>
using not_any_arithmetic_type = disable_if_any<is_arithmetic_decay<Types>...>;

// Check if a type contains or is arithmetic
template <typename T>
struct is_arithmetic_container
    : std::integral_constant<
          bool, std::is_arithmetic<scalar_type_t<std::decay_t<T>>>::value> {};

template <typename T>
using arithmetic_type_container
    = std::enable_if_t<is_arithmetic_container<T>::value>;

template <typename T>
using not_arithmetic_type_container = disable_if<is_arithmetic_container<T>>;

template <typename... Types>
using all_arithmetic_type_container
    = enable_if_all<is_arithmetic_container<Types>...>;

template <typename... Types>
using any_arithmetic_type_container
    = enable_if_any<is_arithmetic_container<Types>...>;

template <typename... Types>
using not_all_arithmetic_type_container
    = disable_if_all<is_arithmetic_container<Types>...>;

template <typename... Types>
using not_any_arithmetic_type_container
    = disable_if_any<is_arithmetic_container<Types>...>;

// Checks whether the type is floating_point
template <typename T>
struct is_fp_decay
    : std::integral_constant<bool,
                             std::is_floating_point<std::decay_t<T>>::value> {};

template <typename T>
using floating_point_type = std::enable_if_t<is_fp_decay<T>::value>;

template <typename T>
using not_floating_point_type = disable_if<is_fp_decay<T>>;

template <typename... Types>
using all_floating_point_type = enable_if_all<is_fp_decay<Types>...>;

template <typename... Types>
using any_floating_point_type = enable_if_any<is_fp_decay<Types>...>;

template <typename... Types>
using not_all_floating_point_type = disable_if_all<is_fp_decay<Types>...>;

template <typename... Types>
using not_any_floating_point_type = disable_if_any<is_fp_decay<Types>...>;

// Check if a type contains or is floating point
template <typename T>
struct is_floating_point_container
    : std::integral_constant<
          bool, std::is_floating_point<scalar_type_t<std::decay_t<T>>>::value> {
};

template <typename T>
using floating_point_type_container
    = std::enable_if_t<is_floating_point_container<T>::value>;

template <typename T>
using not_floating_point_type_container
    = disable_if<is_floating_point_container<T>>;

template <typename... Types>
using all_floating_point_type_container
    = enable_if_all<is_floating_point_container<Types>...>;

template <typename... Types>
using any_floating_point_type_container
    = enable_if_any<is_floating_point_container<Types>...>;

template <typename... Types>
using not_all_floating_point_type_container
    = disable_if_all<is_floating_point_container<Types>...>;

template <typename... Types>
using not_any_floating_point_type_container
    = disable_if_any<is_floating_point_container<Types>...>;

// Check if type is a var

template <typename T>
using var_type = std::enable_if_t<is_var<std::decay<T>>::value>;

template <typename T>
using not_var_type = disable_if<is_var<std::decay<T>>>;

template <typename... Types>
using all_var_type = enable_if_all<is_var<std::decay<Types>>...>;

template <typename... Types>
using any_var_type = enable_if_any<is_var<std::decay<Types>>...>;

template <typename... Types>
using not_all_var_type = disable_if_all<is_var<std::decay<Types>>...>;

template <typename... Types>
using not_any_var_type = disable_if_any<is_var<std::decay<Types>>...>;

// Check if type contains or is a var
template <typename T>
struct is_var_container
    : std::integral_constant<bool,
                             is_var<scalar_type_t<std::decay_t<T>>>::value> {};

template <typename T>
using var_type_container = std::enable_if_t<is_var_container<T>::value>;

template <typename T>
using not_var_type_container = disable_if<is_var_container<T>>;

template <typename... Types>
using all_var_type_container = enable_if_all<is_var_container<Types>...>;

template <typename... Types>
using any_var_type_container = enable_if_any<is_var_container<Types>...>;

template <typename... Types>
using not_all_var_type_container = disable_if_all<is_var_container<Types>...>;

template <typename... Types>
using not_any_var_type_container = disable_if_any<is_var_container<Types>...>;

// Check if type is a fvar
template <typename T>
using fvar_type = std::enable_if_t<is_fvar<std::decay<T>>::value>;

template <typename T>
using not_fvar_type = disable_if<is_fvar<std::decay<T>>>;

template <typename... Types>
using all_fvar_type = enable_if_all<is_fvar<std::decay<Types>>...>;

template <typename... Types>
using any_fvar_type = enable_if_any<is_fvar<std::decay<Types>>...>;

template <typename... Types>
using not_all_fvar_type = disable_if_all<is_fvar<std::decay<Types>>...>;

template <typename... Types>
using not_any_fvar_type = disable_if_any<is_fvar<std::decay<Types>>...>;

// Check if type contains or is a fvar
template <typename T>
struct is_fvar_container
    : std::integral_constant<bool,
                             is_fvar<scalar_type_t<std::decay_t<T>>>::value> {};

template <typename T>
using fvar_type_container = std::enable_if_t<is_fvar_container<T>::value>;

template <typename T>
using not_fvar_type_container = disable_if<is_fvar_container<T>>;

template <typename... Types>
using all_fvar_type_container = enable_if_all<is_fvar_container<Types>...>;

template <typename... Types>
using any_fvar_type_container = enable_if_any<is_fvar_container<Types>...>;

template <typename... Types>
using not_all_fvar_type_container = disable_if_all<is_fvar_container<Types>...>;

template <typename... Types>
using not_any_fvar_type_container = disable_if_any<is_fvar_container<Types>...>;

template <typename T>
struct is_ad_type
    : std::integral_constant<bool, is_fvar<std::decay_t<T>>::value
                                       || is_var<std::decay_t<T>>::value> {};

template <typename T>
using ad_type = std::enable_if_t<is_ad_type<std::decay<T>>::value>;

template <typename T>
using not_ad_type = disable_if<is_ad_type<std::decay<T>>>;

template <typename... Types>
using all_ad_type = enable_if_all<is_ad_type<std::decay<Types>>...>;

template <typename... Types>
using any_ad_type = enable_if_any<is_ad_type<std::decay<Types>>...>;

template <typename... Types>
using not_all_ad_type = disable_if_all<is_ad_type<std::decay<Types>>...>;

template <typename... Types>
using not_any_ad_type = disable_if_any<is_ad_type<std::decay<Types>>...>;

// Check if type contains or is a ad_type
template <typename T>
struct is_ad_type_container
    : std::integral_constant<
          bool, is_ad_type<scalar_type_t<std::decay_t<T>>>::value> {};

template <typename T>
using ad_type_container = std::enable_if_t<is_ad_type_container<T>::value>;

template <typename T>
using not_ad_type_container = disable_if<is_ad_type_container<T>>;

template <typename... Types>
using all_ad_type_container = enable_if_all<is_ad_type_container<Types>...>;

template <typename... Types>
using any_ad_type_container = enable_if_any<is_ad_type_container<Types>...>;

template <typename... Types>
using not_all_ad_type_container
    = disable_if_all<is_ad_type_container<Types>...>;

template <typename... Types>
using not_any_ad_type_container
    = disable_if_any<is_ad_type_container<Types>...>;

// Enables if type is var or arithmetic
template <typename T>
using var_or_arithmetic = std::enable_if_t<is_var_or_arithmetic<T>::value>;

template <typename T>
using not_var_or_arithmetic = disable_if<is_var_or_arithmetic<T>>;

template <typename... Types>
using all_var_or_arithmetic = enable_if_all<is_var_or_arithmetic<Types>...>;

template <typename... Types>
using any_var_or_arithmetic = enable_if_any<is_var_or_arithmetic<Types>...>;

template <typename... Types>
using not_all_var_or_arithmetic
    = disable_if_all<is_var_or_arithmetic<Types>...>;

template <typename... Types>
using not_any_var_or_arithmetic
    = disable_if_any<is_var_or_arithmetic<Types>...>;

// Enables if type is var or arithmetic
template <typename T>
struct is_var_or_arithmetic_container
    : std::integral_constant<bool,
                             is_var_or_arithmetic<scalar_type_t<T>>::value> {};

template <typename T>
using var_or_arithmetic_type_container
    = std::enable_if_t<is_var_or_arithmetic_container<T>::value>;

template <typename T>
using not_var_or_arithmetic_type_container
    = disable_if<is_var_or_arithmetic_container<T>>;

template <typename... Types>
using all_var_or_arithmetic_type_container
    = enable_if_all<is_var_or_arithmetic_container<Types>...>;

template <typename... Types>
using any_var_or_arithmetic_type_container
    = enable_if_any<is_var_or_arithmetic_container<Types>...>;

template <typename... Types>
using not_all_var_or_arithmetic_type_container
    = disable_if_all<is_var_or_arithmetic_container<Types>...>;

template <typename... Types>
using not_any_var_or_arithmetic_type_container
    = disable_if_any<is_var_or_arithmetic_container<Types>...>;

// Checks whether type is arithmetic, var, or fvar
template <typename T>
struct is_stan_scalar
    : std::integral_constant<bool, std::is_arithmetic<std::decay_t<T>>::value
                                       || is_var<std::decay_t<T>>::value
                                       || is_fvar<std::decay_t<T>>::value> {};

template <typename T>
using stan_scalar_type = std::enable_if_t<is_stan_scalar<T>::value>;

template <typename T>
using not_stan_scalar_type = disable_if<is_stan_scalar<T>>;

template <typename... Types>
using all_stan_scalar_type = enable_if_all<is_stan_scalar<Types>...>;

template <typename... Types>
using any_stan_scalar_type = enable_if_any<is_stan_scalar<Types>...>;

template <typename... Types>
using not_all_stan_scalar_type = disable_if_all<is_stan_scalar<Types>...>;

template <typename... Types>
using not_any_stan_scalar_type = disable_if_any<is_stan_scalar<Types>...>;

// Check whether a type contains (or is) arithmetic, var, or fvar
template <typename T>
struct is_stan_scalar_container
    : std::integral_constant<bool, is_stan_scalar<scalar_type_t<T>>::value> {};

template <typename T>
using stan_scalar_type_container
    = std::enable_if_t<is_stan_scalar_container<T>::value>;

template <typename T>
using not_stan_scalar_type_container = disable_if<is_stan_scalar_container<T>>;

template <typename... Types>
using all_stan_scalar_type_container
    = enable_if_all<is_stan_scalar_container<Types>...>;

template <typename... Types>
using any_stan_scalar_type_container
    = enable_if_any<is_stan_scalar_container<Types>...>;

template <typename... Types>
using not_all_stan_scalar_type_container
    = enable_if_all<is_stan_scalar_container<Types>...>;

template <typename... Types>
using not_any_stan_scalar_type_container
    = enable_if_any<is_stan_scalar_container<Types>...>;

// Checks whether type is a scalar as defined by the standard
template <typename T>
struct is_scalar_decay
    : std::integral_constant<bool, std::is_scalar<std::decay_t<T>>::value> {};

// Checks whether decayed type is a vector
template <typename T>
struct is_vector_decay
    : std::integral_constant<bool, is_vector<std::decay_t<T>>::value> {};

template <typename T>
using vector_type = std::enable_if_t<is_vector_decay<T>::value>;

template <typename T>
using not_vector_type = disable_if<is_vector_decay<T>>;

template <typename... Types>
using all_vector_type = enable_if_all<is_vector_decay<Types>...>;

template <typename... Types>
using any_vector_type = enable_if_any<is_vector_decay<Types>...>;

template <typename... Types>
using not_all_vector_type = disable_if_all<is_vector_decay<Types>...>;

template <typename... Types>
using not_any_vector_type = disable_if_any<is_vector_decay<Types>...>;

// Checks whether decayed type is a standard vector
template <typename T>
struct is_std_vector_decay
    : std::integral_constant<bool, is_std_vector<std::decay_t<T>>::value> {};

template <typename T>
using std_vector_type = std::enable_if_t<is_std_vector_decay<T>::value>;

template <typename T>
using not_std_vector_type = disable_if<is_std_vector_decay<T>>;

template <typename... Types>
using all_std_vector_type = enable_if_all<is_std_vector_decay<Types>...>;

template <typename... Types>
using any_std_vector_type = enable_if_any<is_std_vector_decay<Types>...>;

template <typename... Types>
using not_all_std_vector_type = disable_if_all<is_std_vector_decay<Types>...>;

template <typename... Types>
using not_any_std_vector_type = disable_if_any<is_std_vector_decay<Types>...>;

// Enable if function comes from Eigen
template <typename T>
using eigen_type = std::enable_if_t<is_eigen<T>::value>;

template <typename T>
using not_eigen_type = disable_if<is_eigen<T>>;

template <typename... Types>
using all_eigen_type = enable_if_all<is_eigen<Types>...>;

template <typename... Types>
using any_eigen_type = enable_if_any<is_eigen<Types>...>;

template <typename... Types>
using not_all_eigen_type = disable_if_all<is_eigen<Types>...>;

template <typename... Types>
using not_any_eigen_type = disable_if_any<is_eigen<Types>...>;

// Checks if type is Eigen or arithmetic, var, or fvar
template <typename T>
struct is_eigen_or_stan_scalar
    : std::integral_constant<bool,
                             is_eigen<T>::value || is_stan_scalar<T>::value> {};

template <typename T>
using eigen_or_stan_scalar
    = std::enable_if_t<is_eigen_or_stan_scalar<T>::value>;

template <typename T>
using not_eigen_or_stan_scalar_type = disable_if<is_eigen_or_stan_scalar<T>>;

template <typename... Types>
using all_eigen_or_stan_scalar_type
    = enable_if_all<is_eigen_or_stan_scalar<Types>...>;

template <typename... Types>
using any_eigen_or_stan_scalar_type
    = enable_if_any<is_eigen_or_stan_scalar<Types>...>;

template <typename... Types>
using not_all_eigen_or_stan_scalar_type
    = disable_if_all<is_eigen_or_stan_scalar<Types>...>;

template <typename... Types>
using not_any_eigen_or_stan_scalar_type
    = disable_if_any<is_eigen_or_stan_scalar<Types>...>;

// Checks if type is Eigen and scalar type is arithmetic
template <typename T>
struct is_eigen_arithmetic
    : std::integral_constant<bool, is_eigen<T>::value
                                       && is_arithmetic_container<T>::value> {};

template <typename T>
using eigen_arithmetic = std::enable_if_t<is_eigen_arithmetic<T>::value>;

template <typename T>
using not_eigen_arithmetic = disable_if<is_eigen_arithmetic<T>>;

template <typename... Types>
using all_eigen_arithmetic = enable_if_all<is_eigen_arithmetic<Types>...>;

template <typename... Types>
using any_eigen_arithmetic = enable_if_any<is_eigen_arithmetic<Types>...>;

template <typename... Types>
using not_all_eigen_arithmetic = disable_if_all<is_eigen_arithmetic<Types>...>;

template <typename... Types>
using not_any_eigen_arithmetic = disable_if_any<is_eigen_arithmetic<Types>...>;

// Checks if type is Eigen and scalar type is var
template <typename T>
struct is_eigen_var
    : std::integral_constant<bool, is_eigen<T>::value
                                       && is_var_container<T>::value> {};

template <typename T>
using eigen_var = std::enable_if_t<is_eigen_var<T>::value>;

template <typename T>
using not_eigen_var = disable_if<is_eigen_var<T>>;

template <typename... Types>
using all_eigen_var = enable_if_all<is_eigen_var<Types>...>;

template <typename... Types>
using any_eigen_var = enable_if_any<is_eigen_var<Types>...>;

template <typename... Types>
using not_all_eigen_var = disable_if_all<is_eigen_var<Types>...>;

template <typename... Types>
using not_any_eigen_var = disable_if_any<is_eigen_var<Types>...>;

// Checks if type is Eigen and scalar type is fvar
template <typename T>
struct is_eigen_fvar
    : std::integral_constant<bool, is_eigen<T>::value
                                       && is_fvar_container<T>::value> {};

template <typename T>
using eigen_fvar = std::enable_if_t<is_eigen_fvar<T>::value>;

template <typename T>
using not_eigen_fvar = disable_if<is_eigen_fvar<T>>;

template <typename... Types>
using all_eigen_fvar = enable_if_all<is_eigen_fvar<Types>...>;

template <typename... Types>
using any_eigen_fvar = enable_if_any<is_eigen_fvar<Types>...>;

template <typename... Types>
using not_all_eigen_fvar = disable_if_all<is_eigen_fvar<Types>...>;

template <typename... Types>
using not_any_eigen_fvar = disable_if_any<is_eigen_fvar<Types>...>;

// Checks if type is Eigen and scalar type is an autodiff type
template <typename T>
struct is_eigen_ad_type
    : std::integral_constant<bool, is_eigen<T>::value
                                       && is_ad_type_container<T>::value> {};

template <typename T>
using eigen_ad_type = std::enable_if_t<is_eigen_ad_type<T>::value>;

template <typename T>
using not_eigen_ad_type = disable_if<is_eigen_ad_type<T>>;

template <typename... Types>
using all_eigen_ad_type = enable_if_all<is_eigen_ad_type<Types>...>;

template <typename... Types>
using any_eigen_ad_type = enable_if_any<is_eigen_ad_type<Types>...>;

template <typename... Types>
using not_all_eigen_ad_type = disable_if_all<is_eigen_ad_type<Types>...>;

template <typename... Types>
using not_any_eigen_ad_type = disable_if_any<is_eigen_ad_type<Types>...>;

// Enable if for Eigen col vectors
template <typename T>
using eigen_col_vector = std::enable_if_t<is_eigen_col_vector<T>::value>;

template <typename T>
using not_eigen_col_vector = disable_if<is_eigen_col_vector<T>>;

template <typename... Types>
using all_eigen_col_vector = enable_if_all<is_eigen_col_vector<Types>...>;

template <typename... Types>
using any_eigen_col_vector = enable_if_any<is_eigen_col_vector<Types>...>;

template <typename... Types>
using not_all_eigen_col_vector = disable_if_all<is_eigen_col_vector<Types>...>;

template <typename... Types>
using not_any_eigen_col_vector = disable_if_any<is_eigen_col_vector<Types>...>;

// Enable if for Eigen row vectors
template <typename T>
using eigen_row_vector = std::enable_if_t<is_eigen_row_vector<T>::value>;

template <typename T>
using not_eigen_row_vector = disable_if<is_eigen_row_vector<T>>;

template <typename... Types>
using all_eigen_row_vector = enable_if_all<is_eigen_row_vector<Types>...>;

template <typename... Types>
using any_eigen_row_vector = enable_if_any<is_eigen_row_vector<Types>...>;

template <typename... Types>
using not_all_eigen_row_vector = disable_if_all<is_eigen_row_vector<Types>...>;

template <typename... Types>
using not_any_eigen_row_vector = disable_if_any<is_eigen_row_vector<Types>...>;

// Enable if eigen row or column vector
template <typename T>
using eigen_vector = std::enable_if_t<is_eigen_vector<T>::value>;

template <typename T>
using not_eigen_vector = disable_if<is_eigen_vector<T>>;

template <typename... Types>
using all_eigen_vector = enable_if_all<is_eigen_vector<Types>...>;

template <typename... Types>
using any_eigen_vector = enable_if_any<is_eigen_vector<Types>...>;

template <typename... Types>
using not_all_eigen_vector = disable_if_all<is_eigen_vector<Types>...>;

template <typename... Types>
using not_any_eigen_vector = disable_if_any<is_eigen_vector<Types>...>;

// Check whether Eigen types satisfy a dot product
template <typename T1, typename T2>
struct is_dot_product
    : std::integral_constant<bool, is_eigen_row_vector<T1>::value
                                       && is_eigen_col_vector<T2>::value> {};

template <typename T1, typename T2>
using dot_product = std::enable_if_t<is_dot_product<T1, T2>::value>;

template <typename T1, typename T2>
using not_dot_product = disable_if<is_dot_product<T1, T2>>;

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
using eigen_rows_match = std::enable_if_t<is_eigen_rows_match<T1, T2>::value>;

template <typename T1, typename T2>
using not_eigen_rows_match = disable_if<is_eigen_rows_match<T1, T2>>;

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
using eigen_cols_match = std::enable_if_t<is_eigen_cols_match<T1, T2>::value>;

template <typename T1, typename T2>
using not_eigen_cols_match = disable_if<is_eigen_cols_match<T1, T2>>;

// primary template for checking if eigen matrix cols match
template <class T1, class T2>
struct is_eigen_dims_match
    : std::integral_constant<bool, is_eigen_cols_match<T1, T2>::value
                                       && is_eigen_rows_match<T1, T2>::value> {
};

template <typename T1, typename T2>
using eigen_dims_match = std::enable_if_t<is_eigen_dims_match<T1, T2>::value>;

template <typename T1, typename T2>
using not_eigen_dims_match = disable_if<is_eigen_dims_match<T1, T2>>;

}  // namespace stan
#endif
