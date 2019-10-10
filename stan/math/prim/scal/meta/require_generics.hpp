#ifndef STAN_MATH_PRIM_SCAL_META_REQUIRE_GENERICS_HPP
#define STAN_MATH_PRIM_SCAL_META_REQUIRE_GENERICS_HPP

#include <stan/math/prim/scal/meta/bool_constant.hpp>
#include <stan/math/prim/scal/meta/is_eigen.hpp>
#include <stan/math/prim/scal/meta/is_fvar.hpp>
#include <stan/math/prim/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/is_var_or_arithmetic.hpp>
#include <stan/math/prim/scal/meta/is_vector.hpp>
#include <stan/math/prim/scal/meta/is_vector_like.hpp>
#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <stan/math/prim/scal/meta/value_type.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>

#include <type_traits>
#include <string>

namespace stan {

/**
 * \defgroup require_meta require_\* Meta Type Traits
 *
 * The requires type traits are wrappers for turning on and off definitions
 *  during compile time. You can think of these as "legacy concepts" aka I think
 *the closest thing you can get to C++20 concepts without having the requires
 *keyword.
 *
 * Legacy concepts are done through enable_if_t aliases which are called
 *require_*. There are 5 basic types of requires.
 *
 *    requires_*_t: This means the function requires a certain type to turn on.
 *For instance
 *~~~~~{.cpp}
 * // Works for arithmetic types
 * template <typename T, require_arithmetic_t<T>...>
 * auto add(T&& x, T&& y) { return x + y; }
 *~~~~~
 *    require_not_*_t : When we want to disable for a type
 *~~~~~{.cpp}
 * // Works for non-arithmetic types
 * template <typename T, require_not_arithmetic_t<T>...>
 * auto add(T&& x, T&& y) { return x + y; }
 *~~~~~
 *    require_all_*_t : Takes a parameter pack of types to enable if all types
 *satisfy the check.
 *~~~~~{.cpp}
 * // Works if T1 and T2 are arithmetic
 * template <typename T1, typename T2, require_all_arithmetic_t<T1, T2>...>
 * auto add(T1&& x, T2&& y) { return x + y; }
 *~~~~~
 *    require_any_*_t : Takes a parameter pack of types to enable if any of the
 *types satisfy the check.
 *~~~~~{.cpp}
 * // Works if either T1 or T2 enherit from EigenBase
 * template <typename T1, typename T2, require_any_eigen_t<T1, T2>...>
 * auto add(T1&& x, T2&& y) { return x + y; }
 *~~~~~
 *    require_not_any_*_t : Takes a parameter pack of types to enable if any one
 *of the types are not satisfied.
 *~~~~~{.cpp}
 * // Works if either neither T1 or T2 are arithmetic
 * template <typename T1, typename T2, require_not_any_arithmetic_t<T1, T2>...>
 * auto add(T1 x, T2 y) { return x + y; }
 *~~~~~
 *    require_not_all_*_t : Takes a parameter pack of types to enable if all of
 *the types are not satisfied.
 *~~~~~{.cpp}
 * // Works if neither T1 and T2 are arithmetic
 * template <typename T1, typename T2, require_not_all_arithmetic_t<T1, T2>...>
 * auto add(T1 x, T2 y) { return x + y; }
 *~~~~~
 *
 * In addition, std::vector and Eigen types have additional requires to detect
 *if the value_type (the first underlying type) or the scalar_type (the
 *containers underlying scalar type) satisfy a condition to enable a class or
 *function.
 *
 * The container requires have a _vt and _st to symbolize the above. A function
 *that accepts eigen matrices of whose whose value type is floating point types
 *can be defined as
 *
 *~~~~~{.cpp}
 * template <typename Mat1, typename Mat2,
 * require_all_eigen_vt<std::is_floating_point, Mat1, Mat2>...>
 * auto add(Mat1&& A, Mat2&& B) { return A + B;}
 *~~~~~
 *  A function that accepts standard vectors of Eigen vectors whose scalar type
 *is floating point types can be defined as
 *
 *~~~~~{.cpp}
 * template <typename Vec1, typename Vec2,
 *   require_all_std_vector_vt<is_eigen_vector, Vec1, Vec2>...,
 *   require_all_std_vector_st<is_floating_point, Vec1, Vec2>...>
 * auto add(Vec1&& A, Vec2&& B) {
 *   std::vector<decltype<A[0] + B[0]>> return_vec;
 *   std::transform(A.begin(), A.end(), B.begin(), return_vec.begin(),
 *     [](auto&& x, auto&& y) {
 *         return x + y;
 *     });
 *   return return_vec;
 * }
 *~~~~~
 *
 * In general, these methods allow Stan to have more generic types so that the
 * library can forward along Eigen expression and have better move semantics.
 * For instance, the below code will accept any arbitrary eigen expression
 * that, if it's an rvalue, can be forwarded to another function.
 *
 *~~~~~{.cpp}
 *  template <typename Mat1, typename Mat2,
 * require_all_eigen_vt<is_arithmetic, Mat1, Mat2>...>
 * inline auto a_func(Mat1&& m1, Mat2&& m2) {
 *   check_nan(m1);
 *   check_nan(m2);
 *   // If m1 and/or m2 is an rvalue it will be moved over to this function
 *  // instead of copied to an lvalue
 *   auto B = another_func(std::forward<Mat1>(m1), std::forward<Mat2>(m2)); //
 *(3) return B;
 *~~~~~
 */

/** \addtogroup require_meta
 *  @{
 */
/**
 * If condition is true, template is enabled
 */
template <class Check>
using require_t = std::enable_if_t<Check::value>;

/**
 * If condition is false, template is disabled
 */
template <typename Check>
using require_not_t = std::enable_if_t<!Check::value>;

/**
 * If all conditions are true, template is enabled
 * Returns a type void if all conditions are true and otherwise fails.
 */
template <class... Checks>
using require_all_t = std::enable_if_t<math::conjunction<Checks...>::value>;

/**
 * If any condition is true, template is enabled.
 *
 * Returns a type void if any of the conditions are true and otherwise fails.
 */
template <class... Checks>
using require_any_t = std::enable_if_t<math::disjunction<Checks...>::value>;

/**
 * If all conditions are false, template is enabled.
 *
 * Returns a type void if all of the conditions are false.
 */
template <class... Checks>
using require_all_not_t
    = std::enable_if_t<!math::conjunction<Checks...>::value>;

/**
 * If any condition is false, template is enabled.
 *
 * Returns a type void if any of the conditions are false.
 */
template <class... Checks>
using require_any_not_t
    = std::enable_if_t<!math::disjunction<Checks...>::value>;

/**
 * Require both types to be the same
 */
template <typename T, typename S>
using require_same_t
    = require_t<std::is_same<std::decay_t<T>, std::decay_t<S>>>;

/**
 * Require both types to not be the same
 */
template <typename T, typename S>
using require_not_same_t
    = require_not_t<std::is_same<std::decay_t<T>, std::decay_t<S>>>;

/**
 * Require both types to not be the same
 */
template <typename T, typename... Types>
using require_all_same_t
    = require_all_t<std::is_same<std::decay_t<T>, std::decay_t<Types>>...>;

template <typename T, typename... Types>
using require_any_not_same_t
    = require_all_not_t<std::is_same<std::decay_t<T>, std::decay_t<Types>>...>;

/**
 * Require both types to be the same
 */
template <typename T, typename S>
using require_same_st = require_t<std::is_same<scalar_type_t<std::decay_t<T>>,
                                               scalar_type_t<std::decay_t<S>>>>;

/**
 * Require both types to not be the same
 */
template <typename T, typename S>
using require_not_same_st
    = require_not_t<std::is_same<scalar_type_t<std::decay_t<T>>,
                                 scalar_type_t<std::decay_t<S>>>>;

/**
 * Require both types to not be the same
 */
template <typename T, typename... Types>
using require_all_same_st
    = require_all_t<std::is_same<scalar_type_t<std::decay_t<T>>,
                                 scalar_type_t<std::decay_t<Types>>>...>;

/**
 * Require any of the types to not be the same.
 */
template <typename T, typename... Types>
using require_any_not_same_st
    = require_all_not_t<std::is_same<scalar_type_t<std::decay_t<T>>,
                                     scalar_type_t<std::decay_t<Types>>>...>;

/**
 * Require both types to be the same
 */
template <typename T, typename S>
using require_same_vt = require_t<
    std::is_same<value_type_t<std::decay_t<T>>, value_type_t<std::decay_t<S>>>>;

/**
 * Require both types to not be the same
 */
template <typename T, typename S>
using require_not_same_vt = require_not_t<
    std::is_same<value_type_t<std::decay_t<T>>, value_type_t<std::decay_t<S>>>>;

/**
 * Require both types to not be the same
 */
template <typename T, typename... Types>
using require_all_same_vt
    = require_all_t<std::is_same<value_type_t<std::decay_t<T>>,
                                 value_type_t<std::decay_t<Types>>>...>;

template <typename T, typename... Types>
using require_any_not_same_vt
    = require_all_not_t<std::is_same<value_type_t<std::decay_t<T>>,
                                     value_type_t<std::decay_t<Types>>>...>;

/**
 * Require T is convertible to S
 */
template <typename T, typename S>
using require_convertible_t
    = require_t<std::is_convertible<std::decay_t<T>, std::decay_t<S>>>;

/**
 * Require T is not convertible to S
 */
template <typename T, typename S>
using require_not_convertible_t
    = require_not_t<std::is_convertible<std::decay_t<T>, std::decay_t<S>>>;

/**
 * Require T is convertible to all Types
 */
template <typename T, typename... Types>
using require_all_convertible_t = require_all_t<
    std::is_convertible<std::decay_t<T>, std::decay_t<Types>>...>;

/**
 * require T is not convertible to any Types
 */
template <typename T, typename... Types>
using require_any_not_convertible_t = require_all_not_t<
    std::is_convertible<std::decay_t<T>, std::decay_t<Types>>...>;

/**
 * Below are enablers for
 * - Double or Int
 * - Arithmetic
 * - Floating Point
 * - Var
 * - Var or Arithmetic
 * - Fvar
 * - Var or Fvar
 * - Arithmetic Var or Fvar
 */

/**
 * Checks if decayed type is a double or integer
 * @tparam The type to check
 */
template <typename T>
struct is_double_or_int
    : bool_constant<
          math::disjunction<std::is_same<double, std::decay_t<T>>,
                            std::is_same<int, std::decay_t<T>>>::value> {};

template <typename T>
using require_double_or_int_t = require_t<is_double_or_int<std::decay_t<T>>>;

template <typename T>
using require_not_double_or_int_t
    = require_not_t<is_double_or_int<std::decay_t<T>>>;

template <typename... Types>
using require_all_double_or_int_t
    = require_all_t<is_double_or_int<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_double_or_int_t
    = require_any_t<is_double_or_int<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_double_or_int_t
    = require_all_not_t<is_double_or_int<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_double_or_int_t
    = require_any_not_t<is_double_or_int<std::decay_t<Types>>...>;

// Checks for arithmetic types
template <typename T>
using require_arithmetic_t = require_t<std::is_arithmetic<std::decay_t<T>>>;

template <typename T>
using require_not_arithmetic_t
    = require_not_t<std::is_arithmetic<std::decay_t<T>>>;

template <typename... Types>
using require_all_arithmetic_t
    = require_all_t<std::is_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_arithmetic_t
    = require_any_t<std::is_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_arithmetic_t
    = require_all_not_t<std::is_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_arithmetic_t
    = require_any_not_t<std::is_arithmetic<std::decay_t<Types>>...>;

// Checks for floating_point types
template <typename T>
using require_floating_point_t
    = require_t<std::is_floating_point<std::decay_t<T>>>;

template <typename T>
using require_not_floating_point_t
    = require_not_t<std::is_floating_point<std::decay_t<T>>>;

template <typename... Types>
using require_all_floating_point_t
    = require_all_t<std::is_floating_point<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_floating_point_t
    = require_any_t<std::is_floating_point<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_floating_point_t
    = require_all_not_t<std::is_floating_point<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_floating_point_t
    = require_any_not_t<std::is_floating_point<std::decay_t<Types>>...>;

template <typename T>
using require_var_t = require_t<is_var<std::decay_t<T>>>;

template <typename T>
using require_not_var_t = require_not_t<is_var<std::decay_t<T>>>;

template <typename... Types>
using require_all_var_t = require_all_t<is_var<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_var_t = require_any_t<is_var<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_var_t = require_all_not_t<is_var<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_var_t = require_any_not_t<is_var<std::decay_t<Types>>...>;

template <typename T>
using require_var_or_arithmetic_t
    = require_t<is_var_or_arithmetic<std::decay_t<T>>>;

template <typename T>
using require_not_var_or_arithmetic_t
    = require_not_t<is_var_or_arithmetic<std::decay_t<T>>>;

template <typename... Types>
using require_all_var_or_arithmetic_t
    = require_all_t<is_var_or_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_var_or_arithmetic_t
    = require_any_t<is_var_or_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_var_or_arithmetic_t
    = require_all_not_t<is_var_or_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_var_or_arithmetic_t
    = require_any_not_t<is_var_or_arithmetic<std::decay_t<Types>>...>;

template <typename T>
using require_fvar_t = require_t<is_fvar<std::decay_t<T>>>;

template <typename T>
using require_not_fvar_t = require_not_t<is_fvar<std::decay_t<T>>>;

template <typename... Types>
using require_all_fvar_t = require_all_t<is_fvar<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_fvar_t = require_any_t<is_fvar<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_fvar_t
    = require_all_not_t<is_fvar<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_fvar_t
    = require_any_not_t<is_fvar<std::decay_t<Types>>...>;

/**
 * Checks if decayed type is a var or fvar
 * @tparam The type to check
 */
template <typename T>
struct is_var_or_fvar
    : bool_constant<math::disjunction<is_var<std::decay_t<T>>,
                                      is_fvar<std::decay_t<T>>>::value> {};

template <typename T>
using require_var_or_fvar_t = require_t<is_var_or_fvar<T>>;

template <typename T>
using require_not_var_or_fvar_t = require_not_t<is_var_or_fvar<T>>;

template <typename... Types>
using require_all_var_or_fvar_t = require_all_t<is_var_or_fvar<Types>...>;

template <typename... Types>
using require_any_var_or_fvar_t = require_any_t<is_var_or_fvar<Types>...>;

template <typename... Types>
using require_all_not_var_or_fvar_t
    = require_all_not_t<is_var_or_fvar<Types>...>;

template <typename... Types>
using require_any_not_var_or_fvar_t
    = require_any_not_t<is_var_or_fvar<Types>...>;

/**
 * Checks if decayed type is a var, fvar, or arithmetic
 * @tparam The type to check
 */
template <typename T>
struct is_stan_scalar
    : bool_constant<
          math::disjunction<is_var<std::decay_t<T>>, is_fvar<std::decay_t<T>>,
                            std::is_arithmetic<std::decay_t<T>>>::value> {};

template <typename T>
using require_stan_scalar_t = require_t<is_stan_scalar<T>>;

template <typename T>
using require_not_stan_scalar_t = require_not_t<is_stan_scalar<T>>;

template <typename... Types>
using require_all_stan_scalar_t = require_all_t<is_stan_scalar<Types>...>;

template <typename... Types>
using require_any_stan_scalar_t = require_any_t<is_stan_scalar<Types>...>;

template <typename... Types>
using require_all_not_stan_scalar_t
    = require_all_not_t<is_stan_scalar<Types>...>;

template <typename... Types>
using require_any_not_stan_scalar_t
    = require_any_not_t<is_stan_scalar<Types>...>;

/**
 * Requires for containers
 */

template <typename T>
using require_std_vector_t = require_t<is_std_vector<T>>;

template <typename T>
using require_not_std_vector_t = require_not_t<is_std_vector<T>>;

template <typename... Types>
using require_all_std_vector_t = require_all_t<is_std_vector<Types>...>;

template <typename... Types>
using require_any_std_vector_t = require_any_t<is_std_vector<Types>...>;

template <typename... Types>
using require_all_not_std_vector_t = require_all_not_t<is_std_vector<Types>...>;

template <typename... Types>
using require_any_not_std_vector_t = require_any_not_t<is_std_vector<Types>...>;

template <typename T>
using require_vector_t = require_t<is_vector<T>>;

template <typename T>
using require_not_vector_t = require_not_t<is_vector<T>>;

template <typename... Types>
using require_all_vector_t = require_all_t<is_vector<Types>...>;

template <typename... Types>
using require_any_vector_t = require_any_t<is_vector<Types>...>;

template <typename... Types>
using require_all_not_vector_t = require_all_not_t<is_vector<Types>...>;

template <typename... Types>
using require_any_not_vector_t = require_any_not_t<is_vector<Types>...>;

template <typename T>
using require_eigen_t = require_t<is_eigen<T>>;

template <typename T>
using require_not_eigen_t = require_not_t<is_eigen<T>>;

template <typename... Types>
using require_all_eigen_t = require_all_t<is_eigen<Types>...>;

template <typename... Types>
using require_any_eigen_t = require_any_t<is_eigen<Types>...>;

template <typename... Types>
using require_all_not_eigen_t = require_all_not_t<is_eigen<Types>...>;

template <typename... Types>
using require_any_not_eigen_t = require_any_not_t<is_eigen<Types>...>;

template <typename T>
using require_eigen_vector_t = require_t<is_eigen_vector<T>>;

template <typename T>
using require_not_eigen_vector_t = require_not_t<is_eigen_vector<T>>;

template <typename... Types>
using require_all_eigen_vector_t = require_all_t<is_eigen_vector<Types>...>;

template <typename... Types>
using require_any_eigen_vector_t = require_any_t<is_eigen_vector<Types>...>;

template <typename... Types>
using require_all_not_eigen_vector_t
    = require_all_not_t<is_eigen_vector<Types>...>;

template <typename... Types>
using require_any_not_eigen_vector_t
    = require_any_not_t<is_eigen_vector<Types>...>;

template <typename T>
using require_vector_like_t = require_t<is_vector_like<T>>;

template <typename T>
using require_not_vector_like_t = require_not_t<is_vector_like<T>>;

template <typename... Types>
using require_all_vector_like_t = require_all_t<is_vector_like<Types>...>;

template <typename... Types>
using require_any_vector_like_t = require_any_t<is_vector_like<Types>...>;

template <typename... Types>
using require_all_not_vector_like_t
    = require_all_not_t<is_vector_like<Types>...>;

template <typename... Types>
using require_any_not_vector_like_t
    = require_any_not_t<is_vector_like<Types>...>;

/**
 * Used as the base for checking whether a type is a container with
 * an underlying type
 *
 * @tparam ContainerCheck Templated struct or alias that wraps a static constant
 * value called type. Used to check the container satisfies a particular type
 * check.
 * @tparam CheckType Templated struct or alias that wraps a static constant
 * value called type. Used to check the container's underlying type satisfies a
 * particular type check.
 *
 */
template <template <class...> class ContainerCheck,
          template <class...> class TypeCheck, class... Check>
using container_value_type_check_base = bool_constant<
    math::conjunction<ContainerCheck<std::decay_t<Check>>...,
                      TypeCheck<value_type_t<Check>>...>::value>;

/**
 * Check a templated type to see if it and it's inner type pass a condiational
 * test.
 * @tparam ContainerCheck Templated struct or alias that wraps a static constant
 * value called type. Used to check the container satisfies a particular type
 * check. used like template <typename T, require_container_vt<is_std_vector,
 * is_var, T>...>
 */
template <template <class...> class ContainerCheck,
          template <class...> class TypeCheck, class... Check>
using require_container_vt = require_t<
    container_value_type_check_base<ContainerCheck, TypeCheck, Check...>>;

/**
 * std vector
 */
template <template <class...> class TypeCheck, class... Check>
struct is_std_vector_value_check
    : container_value_type_check_base<is_std_vector, TypeCheck, Check...> {};

template <template <class...> class TypeCheck, class... Check>
using require_std_vector_vt
    = require_t<is_std_vector_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_std_vector_vt
    = require_not_t<is_std_vector_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_std_vector_vt
    = require_any_t<is_std_vector_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_std_vector_vt
    = require_any_not_t<is_std_vector_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_std_vector_vt
    = require_all_t<is_std_vector_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_std_vector_vt
    = require_all_not_t<is_std_vector_value_check<TypeCheck, Check>...>;

/**
 * Vectors
 */
template <template <class...> class TypeCheck, class... Check>
struct is_vector_value_check
    : container_value_type_check_base<is_vector, TypeCheck, Check...> {};

template <template <class...> class TypeCheck, class... Check>
using require_vector_vt = require_t<is_vector_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_vector_vt
    = require_not_t<is_vector_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_vector_vt
    = require_any_t<is_vector_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_vector_vt
    = require_any_not_t<is_vector_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_vector_vt
    = require_all_t<is_vector_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_vector_vt
    = require_all_not_t<is_vector_value_check<TypeCheck, Check>...>;

/**
 * Vector Like
 */
template <template <class...> class TypeCheck, class... Check>
struct is_vector_like_value_check
    : container_value_type_check_base<is_vector_like, TypeCheck, Check...> {};

template <template <class...> class TypeCheck, class... Check>
using require_vector_like_vt
    = require_t<is_vector_like_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_vector_like_vt
    = require_not_t<is_vector_like_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_vector_like_vt
    = require_any_t<is_vector_like_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_vector_like_vt
    = require_any_not_t<is_vector_like_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_vector_like_vt
    = require_all_t<is_vector_like_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_vector_like_vt
    = require_all_not_t<is_vector_like_value_check<TypeCheck, Check>...>;

/**
 * Eigen
 */
template <template <class...> class TypeCheck, class... Check>
struct is_eigen_value_check
    : container_value_type_check_base<is_eigen, TypeCheck, Check...> {};

template <template <class...> class TypeCheck, class... Check>
using require_eigen_vt = require_t<is_eigen_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_eigen_vt
    = require_not_t<is_eigen_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_eigen_vt
    = require_any_t<is_eigen_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_eigen_vt
    = require_any_not_t<is_eigen_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_eigen_vt
    = require_all_t<is_eigen_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_eigen_vt
    = require_all_not_t<is_eigen_value_check<TypeCheck, Check>...>;

/**
 * Std vector
 */
template <template <class...> class TypeCheck, class... Check>
struct is_eigen_vector_value_check
    : container_value_type_check_base<is_eigen_vector, TypeCheck, Check...> {};

template <template <class...> class TypeCheck, class... Check>
using require_eigen_vector_vt
    = require_t<is_eigen_vector_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_eigen_vector_vt
    = require_not_t<is_eigen_vector_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_eigen_vector_vt
    = require_any_t<is_eigen_vector_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_eigen_vector_vt
    = require_any_not_t<is_eigen_vector_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_eigen_vector_vt
    = require_all_t<is_eigen_vector_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_eigen_vector_vt
    = require_all_not_t<is_eigen_vector_value_check<TypeCheck, Check>...>;

/**
 * Container Scalar check
 */

/**
 * Used as the base for checking whether a type is a container with
 * an underlying scalar type
 *
 * @tparam ContainerCheck Templated struct or alias that wraps a static constant
 * scalar called type. Used to check the container satisfies a particular type
 * check.
 * @tparam CheckType Templated struct or alias that wraps a static constant
 * scalar called type. Used to check the container's underlying type satisfies a
 * particular type check.
 *
 */
template <template <class...> class ContainerCheck,
          template <class...> class TypeCheck, class... Check>
using container_scalar_type_check_base = bool_constant<
    math::conjunction<ContainerCheck<std::decay_t<Check>>...,
                      TypeCheck<scalar_type_t<Check>>...>::value>;

/**
 * Check a templated type to see if it and it's inner type pass a condiational
 * test.
 * @tparam ContainerCheck Templated struct or alias that wraps a static constant
 * scalar called type. Used to check the container satisfies a particular type
 * check. used like template <typename T, require_container_st<is_std_vector,
 * is_var, T>...>
 */
template <template <class...> class ContainerCheck,
          template <class...> class TypeCheck, class... Check>
using require_container_st = require_t<
    container_scalar_type_check_base<ContainerCheck, TypeCheck, Check...>>;

/**
 * std vector
 */
template <template <class...> class TypeCheck, class... Check>
struct is_std_vector_scalar_check
    : container_scalar_type_check_base<is_std_vector, TypeCheck, Check...> {};

template <template <class...> class TypeCheck, class... Check>
using require_std_vector_st
    = require_t<is_std_vector_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_std_vector_st
    = require_not_t<is_std_vector_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_std_vector_st
    = require_any_t<is_std_vector_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_std_vector_st
    = require_any_not_t<is_std_vector_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_std_vector_st
    = require_all_t<is_std_vector_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_std_vector_st
    = require_all_not_t<is_std_vector_scalar_check<TypeCheck, Check>...>;

/**
 * Vectors
 */
template <template <class...> class TypeCheck, class... Check>
struct is_vector_scalar_check
    : container_scalar_type_check_base<is_vector, TypeCheck, Check...> {};

template <template <class...> class TypeCheck, class... Check>
using require_vector_st
    = require_t<is_vector_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_vector_st
    = require_not_t<is_vector_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_vector_st
    = require_any_t<is_vector_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_vector_st
    = require_any_not_t<is_vector_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_vector_st
    = require_all_t<is_vector_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_vector_st
    = require_all_not_t<is_vector_scalar_check<TypeCheck, Check>...>;

/**
 * Vector Like
 */
template <template <class...> class TypeCheck, class... Check>
struct is_vector_like_scalar_check
    : container_scalar_type_check_base<is_vector_like, TypeCheck, Check...> {};

template <template <class...> class TypeCheck, class... Check>
using require_vector_like_st
    = require_t<is_vector_like_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_vector_like_st
    = require_not_t<is_vector_like_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_vector_like_st
    = require_any_t<is_vector_like_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_vector_like_st
    = require_any_not_t<is_vector_like_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_vector_like_st
    = require_all_t<is_vector_like_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_vector_like_st
    = require_all_not_t<is_vector_like_scalar_check<TypeCheck, Check>...>;

/**
 * Eigen
 */
template <template <class...> class TypeCheck, class... Check>
struct is_eigen_scalar_check
    : container_scalar_type_check_base<is_eigen, TypeCheck, Check...> {};

template <template <class...> class TypeCheck, class... Check>
using require_eigen_st = require_t<is_eigen_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_eigen_st
    = require_not_t<is_eigen_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_eigen_st
    = require_any_t<is_eigen_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_eigen_st
    = require_any_not_t<is_eigen_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_eigen_st
    = require_all_t<is_eigen_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_eigen_st
    = require_all_not_t<is_eigen_scalar_check<TypeCheck, Check>...>;

/**
 * Std vector
 */
template <template <class...> class TypeCheck, class... Check>
struct is_eigen_vector_scalar_check
    : container_scalar_type_check_base<is_eigen_vector, TypeCheck, Check...> {};

template <template <class...> class TypeCheck, class... Check>
using require_eigen_vector_st
    = require_t<is_eigen_vector_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_eigen_vector_st
    = require_not_t<is_eigen_vector_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_eigen_vector_st
    = require_any_t<is_eigen_vector_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_eigen_vector_st
    = require_any_not_t<is_eigen_vector_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_eigen_vector_st
    = require_all_t<is_eigen_vector_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_eigen_vector_st
    = require_all_not_t<is_eigen_vector_scalar_check<TypeCheck, Check>...>;

/** @}*/
}  // namespace stan
#endif
