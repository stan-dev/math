#ifndef STAN_MATH_PRIM_META_REQUIRE_GENERICS_HPP
#define STAN_MATH_PRIM_META_REQUIRE_GENERICS_HPP

#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_container.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_fvar.hpp>
#include <stan/math/prim/meta/is_string_convertible.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/prim/meta/is_var_or_arithmetic.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/is_vector_like.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/conjunction.hpp>
#include <stan/math/prim/meta/disjunction.hpp>

#include <type_traits>
#include <string>

namespace stan {

/** \addtogroup type_trait
 *  @{
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

/**
 * Checks if decayed type is a var or fvar
 * @tparam The type to check
 */
template <typename T>
struct is_autodiff
    : bool_constant<math::disjunction<is_var<std::decay_t<T>>,
                                      is_fvar<std::decay_t<T>>>::value> {};

/**
 * Checks if decayed type is a var, fvar, or arithmetic
 * @tparam The type to check
 */
template <typename T>
struct is_stan_scalar
    : bool_constant<
          math::disjunction<is_var<std::decay_t<T>>, is_fvar<std::decay_t<T>>,
                            std::is_arithmetic<std::decay_t<T>>>::value> {};

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
 * Meta type trait to check if type is a standard vector and value type passes
 * additional type trait check
 */
template <template <class...> class TypeCheck, class... Check>
struct is_std_vector_value_check
    : container_value_type_check_base<is_std_vector, TypeCheck, Check...> {};

/**
 * Meta type trait to check if type is a standard vector or eigen vector and
 * value type passes additional type trait check
 */
template <template <class...> class TypeCheck, class... Check>
struct is_vector_value_check
    : container_value_type_check_base<is_vector, TypeCheck, Check...> {};

/**
 * Meta type trait to check if type is a standard vector, eigen vector, or array
 * and value type passes additional type trait check
 */
template <template <class...> class TypeCheck, class... Check>
struct is_vector_like_value_check
    : container_value_type_check_base<is_vector_like, TypeCheck, Check...> {};

/**
 * Meta type trait to check if type inherits from EigenBase.
 */
template <template <class...> class TypeCheck, class... Check>
struct is_eigen_value_check
    : container_value_type_check_base<is_eigen, TypeCheck, Check...> {};

/**
 * Meta type trait to check if type inherits from EigenBase and has row or col
 * values at compile time of 1.
 */
template <template <class...> class TypeCheck, class... Check>
struct is_eigen_vector_value_check
    : container_value_type_check_base<is_eigen_vector, TypeCheck, Check...> {};

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
 * Meta type trait to check if type is a standard vector and scalar type passes
 * additional type trait check
 */
template <template <class...> class TypeCheck, class... Check>
struct is_std_vector_scalar_check
    : container_scalar_type_check_base<is_std_vector, TypeCheck, Check...> {};

/**
 * Meta type trait to check if type is a standard vector or eigen vector and
 * scalar type passes additional type trait check
 */
template <template <class...> class TypeCheck, class... Check>
struct is_vector_scalar_check
    : container_scalar_type_check_base<is_vector, TypeCheck, Check...> {};

/**
 * Meta type trait to check if type is a standard vector, eigen vector, or array
 * and scalar type passes additional type trait check
 */
template <template <class...> class TypeCheck, class... Check>
struct is_vector_like_scalar_check
    : container_scalar_type_check_base<is_vector_like, TypeCheck, Check...> {};

/**
 * Meta type trait to check if type inherits from EigenBase and scalar type pass
 * additional type trait check
 */
template <template <class...> class TypeCheck, class... Check>
struct is_eigen_scalar_check
    : container_scalar_type_check_base<is_eigen, TypeCheck, Check...> {};

/**
 * Meta type trait to check if type inherits from EigenBase, has 1 row or column
 * compile time size, and scalar type passes additional type trait check
 */
template <template <class...> class TypeCheck, class... Check>
struct is_eigen_vector_scalar_check
    : container_scalar_type_check_base<is_eigen_vector, TypeCheck, Check...> {};
/** @}*/

/** \addtogroup require_base_types
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
    = std::enable_if_t<!math::disjunction<Checks...>::value>;

/**
 * If any condition is false, template is enabled.
 *
 * Returns a type void if any of the conditions are false.
 */
template <class... Checks>
using require_any_not_t
    = std::enable_if_t<!math::conjunction<Checks...>::value>;

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

/**
 * Require the first type to differ from at least one of the rest
 */
template <typename T, typename... Types>
using require_any_not_same_t
    = require_any_not_t<std::is_same<std::decay_t<T>, std::decay_t<Types>>...>;

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
    = require_any_not_t<std::is_same<scalar_type_t<std::decay_t<T>>,
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
    = require_any_not_t<std::is_same<value_type_t<std::decay_t<T>>,
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
using require_any_not_convertible_t = require_any_not_t<
    std::is_convertible<std::decay_t<T>, std::decay_t<Types>>...>;

/**
 * Checks if type is implicitly convertible to std::string
 */
template <typename T>
using require_string_convertible_t
    = require_t<is_string_convertible<std::decay_t<T>>>;

template <typename T>
using require_not_string_convertible_t
    = require_not_t<is_string_convertible<std::decay_t<T>>>;

template <typename... Types>
using require_all_string_convertible_t
    = require_all_t<is_string_convertible<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_string_convertible_t
    = require_any_t<is_string_convertible<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_string_convertible_t
    = require_all_not_t<is_string_convertible<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_string_convertible_t
    = require_any_not_t<is_string_convertible<std::decay_t<Types>>...>;

/**
 * Below are enablers for
 * - Double or Int
 * - Arithmetic
 * - Floating Point
 * - Index
 * - Var
 * - Var or Arithmetic
 * - Fvar
 * - Var or Fvar
 * - Arithmetic Var or Fvar
 */

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

// Checks if type is something we would use for index
template <typename T>
using require_index_t = require_t<std::is_integral<std::decay_t<T>>>;

template <typename T>
using require_not_index_t = require_not_t<std::is_integral<std::decay_t<T>>>;

template <typename... Types>
using require_all_index_t
    = require_all_t<std::is_integral<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_index_t
    = require_any_t<std::is_integral<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_index_t
    = require_all_not_t<std::is_integral<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_index_t
    = require_any_not_t<std::is_integral<std::decay_t<Types>>...>;

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

template <typename T>
using require_autodiff_t = require_t<is_autodiff<T>>;

template <typename T>
using require_not_autodiff_t = require_not_t<is_autodiff<T>>;

template <typename... Types>
using require_all_autodiff_t = require_all_t<is_autodiff<Types>...>;

template <typename... Types>
using require_any_autodiff_t = require_any_t<is_autodiff<Types>...>;

template <typename... Types>
using require_all_not_autodiff_t = require_all_not_t<is_autodiff<Types>...>;

template <typename... Types>
using require_any_not_autodiff_t = require_any_not_t<is_autodiff<Types>...>;

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
/** @}*/

/**
 * Requires for containers
 */

/** \addtogroup require_container_types
 *  @{
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

// Checks for container types
template <typename T>
using require_container_t = require_t<is_container<std::decay_t<T>>>;

template <typename T>
using require_not_container_t = require_not_t<is_container<std::decay_t<T>>>;

template <typename... Types>
using require_all_container_t
    = require_all_t<is_container<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_container_t
    = require_any_t<is_container<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_container_t
    = require_all_not_t<is_container<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_container_t
    = require_any_not_t<is_container<std::decay_t<Types>>...>;

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
 * Check a templated type to see if it and its inner type pass a conditional
 * test.
 * @tparam ContainerCheck Templated struct or alias that wraps a static constant
 * scalar called type. Used to check the container satisfies a particular type
 * check. Used like template <typename T, require_container_st<is_std_vector,
 * is_var, T>...>
 */
template <template <class...> class ContainerCheck,
          template <class...> class TypeCheck, class... Check>
using require_container_st = require_t<
    container_scalar_type_check_base<ContainerCheck, TypeCheck, Check...>>;

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

template <typename Row, typename Col>
using require_eigen_row_and_col_t = require_t<
    math::conjunction<is_eigen_row_vector<Row>, is_eigen_col_vector<Col>>>;

template <typename Row, typename Col>
using require_not_eigen_row_and_col_t = require_not_t<
    math::conjunction<is_eigen_row_vector<Row>, is_eigen_col_vector<Col>>>;

/** @}*/
}  // namespace stan
#endif
