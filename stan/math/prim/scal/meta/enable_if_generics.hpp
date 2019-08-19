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
using DisableIf = typename std::enable_if_t<!Check::value>;

/**
 * Template metaprogram for deducing if a method should be remove a
 * a function from overload resolution.
 *
 * Returns a type void if all conditions are true and otherwise fails.
 */
template <class... Checks>
using EnableIfAll = std::enable_if_t<math::conjunction<Checks...>::value>;

/**
 * Template metaprogram for deducing if a method should be remove a
 * a function from overload resolution.
 *
 * Returns a type void if any of the conditions are true and otherwise fails.
 */
template <class... Checks>
using EnableIfAny = std::enable_if_t<math::disjunction<Checks...>::value>;

/**
 * Template metaprogram for deducing if a method should be remove a
 * a function from overload resolution.
 *
 * Returns a type void if all of the conditions are false.
 */
template <class... Checks>
using DisableIfAll = std::enable_if_t<!math::disjunction<Checks...>::value>;

/**
 * Template metaprogram for deducing if a method should be remove a
 * a function from overload resolution.
 *
 * Returns a type void if any of the conditions are false.
 */
template <class... Checks>
using DisableIfAny = std::enable_if_t<!math::conjunction<Checks...>::value>;

// Check whether the decayed types are the same
template <typename T, typename S>
struct is_same_decay
    : std::integral_constant<
          bool, std::is_same<std::decay_t<T>, std::decay_t<S>>::value> {};

template <typename T, typename S>
using SameType = std::enable_if_t<is_same_decay<T, S>::value>;

template <typename T, typename S>
using NotSameType = DisableIf<is_same_decay<T, S>>;

template <typename T, typename... Types>
using AllSameType = EnableIfAll<is_same_decay<T, Types>...>;

template <typename T, typename... Types>
using NotAllSameType = DisableIfAll<is_same_decay<T, Types>...>;

// Checks decayed (non-const/ref'd) type is arithmetic
template <typename T>
struct is_arithmetic_decay
    : std::integral_constant<bool, std::is_arithmetic<std::decay_t<T>>::value> {
};

template <typename T>
using ArithmeticType = std::enable_if_t<is_arithmetic_decay<T>::value>;

template <typename T>
using NotArithmeticType = DisableIf<is_arithmetic_decay<T>>;

template <typename... Types>
using AllArithmeticType = EnableIfAll<is_arithmetic_decay<Types>...>;

template <typename... Types>
using AnyArithmeticType = EnableIfAny<is_arithmetic_decay<Types>...>;

template <typename... Types>
using NotAllArithmeticType = DisableIfAll<is_arithmetic_decay<Types>...>;

template <typename... Types>
using NotAnyArithmeticType = DisableIfAny<is_arithmetic_decay<Types>...>;

// Check if a type contains or is arithmetic
template <typename T>
struct is_arithmetic_container
    : std::integral_constant<
          bool, std::is_arithmetic<scalar_type_t<std::decay_t<T>>>::value> {};

template <typename T>
using ArithmeticTypeContainer
    = std::enable_if_t<is_arithmetic_container<T>::value>;

template <typename T>
using NotArithmeticTypeContainer = DisableIf<is_arithmetic_container<T>>;

template <typename... Types>
using AllArithmeticTypeContainer
    = EnableIfAll<is_arithmetic_container<Types>...>;

template <typename... Types>
using AnyArithmeticTypeContainer
    = EnableIfAny<is_arithmetic_container<Types>...>;

template <typename... Types>
using NotAllArithmeticTypeContainer
    = DisableIfAll<is_arithmetic_container<Types>...>;

template <typename... Types>
using NotAnyArithmeticTypeContainer
    = DisableIfAny<is_arithmetic_container<Types>...>;

// Checks whether the type is floating_point
template <typename T>
struct is_fp_decay
    : std::integral_constant<bool,
                             std::is_floating_point<std::decay_t<T>>::value> {};

template <typename T>
using FloatingPoint = std::enable_if_t<is_fp_decay<T>::value>;

template <typename T>
using NotFloatingPoint = DisableIf<is_fp_decay<T>>;

template <typename... Types>
using AllFloatingPoint = EnableIfAll<is_fp_decay<Types>...>;

template <typename... Types>
using AnyFloatingPoint = EnableIfAny<is_fp_decay<Types>...>;

template <typename... Types>
using NotAllFloatingPoint = DisableIfAll<is_fp_decay<Types>...>;

template <typename... Types>
using NotAnyFloatingPoint = DisableIfAny<is_fp_decay<Types>...>;

// Check if a type contains or is floating point
template <typename T>
struct is_floating_point_container
    : std::integral_constant<
          bool, std::is_floating_point<scalar_type_t<std::decay_t<T>>>::value> {
};

template <typename T>
using FloatingPointTypeContainer
    = std::enable_if_t<is_floating_point_container<T>::value>;

template <typename T>
using NotFloatingPointTypeContainer = DisableIf<is_floating_point_container<T>>;

template <typename... Types>
using AllFloatingPointTypeContainer
    = EnableIfAll<is_floating_point_container<Types>...>;

template <typename... Types>
using AnyFloatingPointTypeContainer
    = EnableIfAny<is_floating_point_container<Types>...>;

template <typename... Types>
using NotAllFloatingPointTypeContainer
    = DisableIfAll<is_floating_point_container<Types>...>;

template <typename... Types>
using NotAnyFloatingPointTypeContainer
    = DisableIfAny<is_floating_point_container<Types>...>;

// Check if type is a var

template <typename T>
using VarType = std::enable_if_t<is_var<std::decay<T>>::value>;

template <typename T>
using NotVarType = DisableIf<is_var<std::decay<T>>>;

template <typename... Types>
using AllVarType = EnableIfAll<is_var<std::decay<Types>>...>;

template <typename... Types>
using AnyVarType = EnableIfAny<is_var<std::decay<Types>>...>;

template <typename... Types>
using NotAllVarType = DisableIfAll<is_var<std::decay<Types>>...>;

template <typename... Types>
using NotAnyVarType = DisableIfAny<is_var<std::decay<Types>>...>;

// Check if type contains or is a var
template <typename T>
struct is_var_container
    : std::integral_constant<bool,
                             is_var<scalar_type_t<std::decay_t<T>>>::value> {};

template <typename T>
using VarTypeContainer = std::enable_if_t<is_var_container<T>::value>;

template <typename T>
using NotVarTypeContainer = DisableIf<is_var_container<T>>;

template <typename... Types>
using AllVarTypeContainer = EnableIfAll<is_var_container<Types>...>;

template <typename... Types>
using AnyVarTypeContainer = EnableIfAny<is_var_container<Types>...>;

template <typename... Types>
using NotAllVarTypeContainer = DisableIfAll<is_var_container<Types>...>;

template <typename... Types>
using NotAnyVarTypeContainer = DisableIfAny<is_var_container<Types>...>;

// Check if type is a fvar
template <typename T>
using FvarType = std::enable_if_t<is_fvar<std::decay<T>>::value>;

template <typename T>
using NotFvarType = DisableIf<is_fvar<std::decay<T>>>;

template <typename... Types>
using AllFvarType = EnableIfAll<is_fvar<std::decay<Types>>...>;

template <typename... Types>
using AnyFvarType = EnableIfAny<is_fvar<std::decay<Types>>...>;

template <typename... Types>
using NotAllFvarType = DisableIfAll<is_fvar<std::decay<Types>>...>;

template <typename... Types>
using NotAnyFvarType = DisableIfAny<is_fvar<std::decay<Types>>...>;

// Check if type contains or is a fvar
template <typename T>
struct is_fvar_container
    : std::integral_constant<bool,
                             is_fvar<scalar_type_t<std::decay_t<T>>>::value> {};

template <typename T>
using fvarFvarTypeContainer = std::enable_if_t<is_fvar_container<T>::value>;

template <typename T>
using NotFvarTypeContainer = DisableIf<is_fvar_container<T>>;

template <typename... Types>
using AllFvarTypeContainer = EnableIfAll<is_fvar_container<Types>...>;

template <typename... Types>
using AnyFvarTypeContainer = EnableIfAny<is_fvar_container<Types>...>;

template <typename... Types>
using NotAllFvarTypeContainer = DisableIfAll<is_fvar_container<Types>...>;

template <typename... Types>
using NotAnyFvarTypeContainer = DisableIfAny<is_fvar_container<Types>...>;

template <typename T>
struct is_ad_type
    : std::integral_constant<bool, is_fvar<std::decay_t<T>>::value
                                       || is_var<std::decay_t<T>>::value> {};

template <typename T>
using AdType = std::enable_if_t<is_ad_type<std::decay<T>>::value>;

template <typename T>
using NotAdType = DisableIf<is_ad_type<std::decay<T>>>;

template <typename... Types>
using AllAdType = EnableIfAll<is_ad_type<std::decay<Types>>...>;

template <typename... Types>
using AnyAdType = EnableIfAny<is_ad_type<std::decay<Types>>...>;

template <typename... Types>
using NotAllAdType = DisableIfAll<is_ad_type<std::decay<Types>>...>;

template <typename... Types>
using NotAnyAdType = DisableIfAny<is_ad_type<std::decay<Types>>...>;

// Check if type contains or is a ad_type
template <typename T>
struct is_ad_type_container
    : std::integral_constant<
          bool, is_ad_type<scalar_type_t<std::decay_t<T>>>::value> {};

template <typename T>
using AdTypeContainer = std::enable_if_t<is_ad_type_container<T>::value>;

template <typename T>
using NotAdTypeContainer = DisableIf<is_ad_type_container<T>>;

template <typename... Types>
using AllAdTypeContainer = EnableIfAll<is_ad_type_container<Types>...>;

template <typename... Types>
using AnyAdTypeContainer = EnableIfAny<is_ad_type_container<Types>...>;

template <typename... Types>
using NotAllAdTypeContainer = DisableIfAll<is_ad_type_container<Types>...>;

template <typename... Types>
using NotAnyAdTypeContainer = DisableIfAny<is_ad_type_container<Types>...>;

// Enables if type is var or arithmetic
template <typename T>
using VarOrArithmetic = std::enable_if_t<is_var_or_arithmetic<T>::value>;

template <typename T>
using NotVarOrArithmetic = DisableIf<is_var_or_arithmetic<T>>;

template <typename... Types>
using AllVarOrArithmetic = EnableIfAll<is_var_or_arithmetic<Types>...>;

template <typename... Types>
using AnyVarOrArithmetic = EnableIfAny<is_var_or_arithmetic<Types>...>;

template <typename... Types>
using NotAllVarOrArithmetic = DisableIfAll<is_var_or_arithmetic<Types>...>;

template <typename... Types>
using NotAnyVarOrArithmetic = DisableIfAny<is_var_or_arithmetic<Types>...>;

// Enables if type is var or arithmetic
template <typename T>
struct is_var_or_arithmetic_container
    : std::integral_constant<bool,
                             is_var_or_arithmetic<scalar_type_t<T>>::value> {};

template <typename T>
using VarOrArithmeticTypeContainer
    = std::enable_if_t<is_var_or_arithmetic_container<T>::value>;

template <typename T>
using NotVarOrArithmeticTypeContainer
    = DisableIf<is_var_or_arithmetic_container<T>>;

template <typename... Types>
using AllVarOrArithmeticTypeContainer
    = EnableIfAll<is_var_or_arithmetic_container<Types>...>;

template <typename... Types>
using AnyVarOrArithmeticTypeContainer
    = EnableIfAny<is_var_or_arithmetic_container<Types>...>;

template <typename... Types>
using NotAllVarOrArithmeticTypeContainer
    = DisableIfAll<is_var_or_arithmetic_container<Types>...>;

template <typename... Types>
using NotAnyVarOrArithmeticTypeContainer
    = DisableIfAny<is_var_or_arithmetic_container<Types>...>;

// Checks whether type is arithmetic, var, or fvar
template <typename T>
struct is_stan_scalar
    : std::integral_constant<bool, std::is_arithmetic<std::decay_t<T>>::value
                                       || is_var<std::decay_t<T>>::value
                                       || is_fvar<std::decay_t<T>>::value> {};

template <typename T>
using StanScalarType = std::enable_if_t<is_stan_scalar<T>::value>;

template <typename T>
using NotStanScalarType = DisableIf<is_stan_scalar<T>>;

template <typename... Types>
using AllStanScalarType = EnableIfAll<is_stan_scalar<Types>...>;

template <typename... Types>
using AnyStanScalarType = EnableIfAny<is_stan_scalar<Types>...>;

template <typename... Types>
using NotAllStanScalarType = DisableIfAll<is_stan_scalar<Types>...>;

template <typename... Types>
using NotAnyStanScalarType = DisableIfAny<is_stan_scalar<Types>...>;

// Check whether a type contains (or is) arithmetic, var, or fvar
template <typename T>
struct is_stan_scalar_container
    : std::integral_constant<bool, is_stan_scalar<scalar_type_t<T>>::value> {};

template <typename T>
using StanScalarTypeContainer
    = std::enable_if_t<is_stan_scalar_container<T>::value>;

template <typename T>
using NotStanScalarTypeContainer = DisableIf<is_stan_scalar_container<T>>;

template <typename... Types>
using AllStanScalarTypeContainer
    = EnableIfAll<is_stan_scalar_container<Types>...>;

template <typename... Types>
using AnyStanScalarTypeContainer
    = EnableIfAny<is_stan_scalar_container<Types>...>;

template <typename... Types>
using NotAllStanScalarTypeContainer
    = EnableIfAll<is_stan_scalar_container<Types>...>;

template <typename... Types>
using NotAnyStanScalarTypeContainer
    = EnableIfAny<is_stan_scalar_container<Types>...>;

// Checks whether type is a scalar as defined by the standard
template <typename T>
struct is_scalar_decay
    : std::integral_constant<bool, std::is_scalar<std::decay_t<T>>::value> {};

template <typename T>
using ScalarType = std::enable_if_t<is_scalar_decay<T>::value>;

template <typename T>
using NotScalarType = DisableIf<is_scalar_decay<T>>;

template <typename... Types>
using AllScalar = EnableIfAll<is_scalar_decay<Types>...>;

template <typename... Types>
using AnyScalar = EnableIfAny<is_scalar_decay<Types>...>;

template <typename... Types>
using NotAllScalar = DisableIfAll<is_scalar_decay<Types>...>;

template <typename... Types>
using NotAnyScalar = DisableIfAny<is_scalar_decay<Types>...>;

// Checks whether decayed type is a vector
template <typename T>
struct is_vector_decay
    : std::integral_constant<bool, is_vector<std::decay_t<T>>::value> {};

template <typename T>
using VectorType = std::enable_if_t<is_vector_decay<T>::value>;

template <typename T>
using NotVectorType = DisableIf<is_vector_decay<T>>;

template <typename... Types>
using AllVectorType = EnableIfAll<is_vector_decay<Types>...>;

template <typename... Types>
using AnyVectorType = EnableIfAny<is_vector_decay<Types>...>;

template <typename... Types>
using NotAllVectorType = DisableIfAll<is_vector_decay<Types>...>;

template <typename... Types>
using NotAnyVectorType = DisableIfAny<is_vector_decay<Types>...>;

// Checks whether decayed type is a standard vector
template <typename T>
struct is_std_vector_decay
    : std::integral_constant<bool, is_std_vector<std::decay_t<T>>::value> {};

template <typename T>
using StdVectorType = std::enable_if_t<is_std_vector_decay<T>::value>;

template <typename T>
using NotStdVectorType = DisableIf<is_std_vector_decay<T>>;

template <typename... Types>
using AllStdVectorType = EnableIfAll<is_std_vector_decay<Types>...>;

template <typename... Types>
using AnyStdVectorType = EnableIfAny<is_std_vector_decay<Types>...>;

template <typename... Types>
using NotAllStdVectorType = DisableIfAll<is_std_vector_decay<Types>...>;

template <typename... Types>
using NotAnyStdVectorType = DisableIfAny<is_std_vector_decay<Types>...>;

// Enable if function comes from Eigen
template <typename T>
using EigenType = std::enable_if_t<is_eigen<T>::value>;

template <typename T>
using NotEigenType = DisableIf<is_eigen<T>>;

template <typename... Types>
using AllEigenType = EnableIfAll<is_eigen<Types>...>;

template <typename... Types>
using AnyEigenType = EnableIfAny<is_eigen<Types>...>;

template <typename... Types>
using NotAllEigenType = DisableIfAll<is_eigen<Types>...>;

template <typename... Types>
using NotAnyEigenType = DisableIfAny<is_eigen<Types>...>;

// Checks if type is Eigen or arithmetic, var, or fvar
template <typename T>
struct is_eigen_or_stan_scalar
    : std::integral_constant<bool,
                             is_eigen<T>::value || is_stan_scalar<T>::value> {};

template <typename T>
using Eigen_or_stan_scalar
    = std::enable_if_t<is_eigen_or_stan_scalar<T>::value>;

template <typename T>
using NotEigenOrStanScalarType = DisableIf<is_eigen_or_stan_scalar<T>>;

template <typename... Types>
using AllEigenOrStanScalarType = EnableIfAll<is_eigen_or_stan_scalar<Types>...>;

template <typename... Types>
using AnyEigenOrStanScalarType = EnableIfAny<is_eigen_or_stan_scalar<Types>...>;

template <typename... Types>
using NotAllEigenOrStanScalarType
    = DisableIfAll<is_eigen_or_stan_scalar<Types>...>;

template <typename... Types>
using NotAnyEigenOrStanScalarType
    = DisableIfAny<is_eigen_or_stan_scalar<Types>...>;

// Checks if type is Eigen and scalar type is arithmetic
template <typename T>
struct is_eigen_arithmetic
    : std::integral_constant<bool, is_eigen<T>::value
                                       && is_arithmetic_container<T>::value> {};

template <typename T>
using EigenArithmetic = std::enable_if_t<is_eigen_arithmetic<T>::value>;

template <typename T>
using NotEigenArithmetic = DisableIf<is_eigen_arithmetic<T>>;

template <typename... Types>
using AllEigenArithmetic = EnableIfAll<is_eigen_arithmetic<Types>...>;

template <typename... Types>
using AnyEigenArithmetic = EnableIfAny<is_eigen_arithmetic<Types>...>;

template <typename... Types>
using NotAllEigenArithmetic = DisableIfAll<is_eigen_arithmetic<Types>...>;

template <typename... Types>
using NotAnyEigenArithmetic = DisableIfAny<is_eigen_arithmetic<Types>...>;

// Checks if type is Eigen and scalar type is var
template <typename T>
struct is_eigen_var
    : std::integral_constant<bool, is_eigen<T>::value
                                       && is_var_container<T>::value> {};

template <typename T>
using EigenVar = std::enable_if_t<is_eigen_var<T>::value>;

template <typename T>
using NotEigenVar = DisableIf<is_eigen_var<T>>;

template <typename... Types>
using AllEigenVar = EnableIfAll<is_eigen_var<Types>...>;

template <typename... Types>
using AnyEigenVar = EnableIfAny<is_eigen_var<Types>...>;

template <typename... Types>
using NotAllEigenVar = DisableIfAll<is_eigen_var<Types>...>;

template <typename... Types>
using NotAnyEigenVar = DisableIfAny<is_eigen_var<Types>...>;

// Checks if type is Eigen and scalar type is fvar
template <typename T>
struct is_eigen_fvar
    : std::integral_constant<bool, is_eigen<T>::value
                                       && is_fvar_container<T>::value> {};

template <typename T>
using EigenFvar = std::enable_if_t<is_eigen_fvar<T>::value>;

template <typename T>
using NotEigenFvar = DisableIf<is_eigen_fvar<T>>;

template <typename... Types>
using AllEigenFvar = EnableIfAll<is_eigen_fvar<Types>...>;

template <typename... Types>
using AnyEigenFvar = EnableIfAny<is_eigen_fvar<Types>...>;

template <typename... Types>
using NotAllEigenFvar = DisableIfAll<is_eigen_fvar<Types>...>;

template <typename... Types>
using NotAnyEigenFvar = DisableIfAny<is_eigen_fvar<Types>...>;

// Checks if type is Eigen and scalar type is an autodiff type
template <typename T>
struct is_eigen_ad_type
    : std::integral_constant<bool, is_eigen<T>::value
                                       && is_ad_type_container<T>::value> {};

template <typename T>
using EigenAdType = std::enable_if_t<is_eigen_ad_type<T>::value>;

template <typename T>
using NotEigenAdType = DisableIf<is_eigen_ad_type<T>>;

template <typename... Types>
using AllEigenAdType = EnableIfAll<is_eigen_ad_type<Types>...>;

template <typename... Types>
using AnyEigenAdType = EnableIfAny<is_eigen_ad_type<Types>...>;

template <typename... Types>
using NotAllEigenAdType = DisableIfAll<is_eigen_ad_type<Types>...>;

template <typename... Types>
using NotAnyEigenAdType = DisableIfAny<is_eigen_ad_type<Types>...>;

// Enable if for Eigen col vectors
template <typename T>
using Eigen_col_vector = std::enable_if_t<is_eigen_col_vector<T>::value>;

template <typename T>
using NotEigenColVector = DisableIf<is_eigen_col_vector<T>>;

template <typename... Types>
using AllEigenColVector = EnableIfAll<is_eigen_col_vector<Types>...>;

template <typename... Types>
using AnyEigenColVector = EnableIfAny<is_eigen_col_vector<Types>...>;

template <typename... Types>
using NotAllEigenColVector = DisableIfAll<is_eigen_col_vector<Types>...>;

template <typename... Types>
using NotAnyEigenColVector = DisableIfAny<is_eigen_col_vector<Types>...>;

// Enable if for Eigen row vectors
template <typename T>
using EigenRowVector = std::enable_if_t<is_eigen_row_vector<T>::value>;

template <typename T>
using NotEigenRowVector = DisableIf<is_eigen_row_vector<T>>;

template <typename... Types>
using AllEigenRowVector = EnableIfAll<is_eigen_row_vector<Types>...>;

template <typename... Types>
using AnyEigenRowVector = EnableIfAny<is_eigen_row_vector<Types>...>;

template <typename... Types>
using NotAllEigenRowVector = DisableIfAll<is_eigen_row_vector<Types>...>;

template <typename... Types>
using NotAnyEigenRowVector = DisableIfAny<is_eigen_row_vector<Types>...>;

// Enable if eigen row or column vector
template <typename T>
using EigenVector = std::enable_if_t<is_eigen_vector<T>::value>;

template <typename T>
using NotEigenVector = DisableIf<is_eigen_vector<T>>;

template <typename... Types>
using AllEigenVector = EnableIfAll<is_eigen_vector<Types>...>;

template <typename... Types>
using AnyEigenVector = EnableIfAny<is_eigen_vector<Types>...>;

template <typename... Types>
using NotAllEigenVector = DisableIfAll<is_eigen_vector<Types>...>;

template <typename... Types>
using NotAnyEigenVector = DisableIfAny<is_eigen_vector<Types>...>;

// Check whether Eigen types satisfy a dot product
template <typename T1, typename T2>
struct is_dot_product
    : std::integral_constant<bool, is_eigen_row_vector<T1>::value
                                       && is_eigen_col_vector<T2>::value> {};

template <typename T1, typename T2>
using DotProduct = std::enable_if_t<is_dot_product<T1, T2>::value>;

template <typename T1, typename T2>
using NotDotProduct = DisableIf<is_dot_product<T1, T2>>;

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
using EigenRowsMatch = std::enable_if_t<is_eigen_rows_match<T1, T2>::value>;

template <typename T1, typename T2>
using NotEigenRowsMatch = DisableIf<is_eigen_rows_match<T1, T2>>;

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
using EigenColsMatch = std::enable_if_t<is_eigen_cols_match<T1, T2>::value>;

template <typename T1, typename T2>
using NotEigenColsMatch = DisableIf<is_eigen_cols_match<T1, T2>>;

// primary template for checking if eigen matrix cols match
template <class T1, class T2>
struct is_eigen_dims_match
    : std::integral_constant<bool, is_eigen_cols_match<T1, T2>::value
                                       && is_eigen_rows_match<T1, T2>::value> {
};

template <typename T1, typename T2>
using EigenDimsMatch = std::enable_if_t<is_eigen_dims_match<T1, T2>::value>;

template <typename T1, typename T2>
using NotEigenDimsMatch = DisableIf<is_eigen_dims_match<T1, T2>>;

}  // namespace stan
#endif
