#ifndef STAN_MATH_PRIM_SCAL_META_REQUIRE_GENERICS_HPP
#define STAN_MATH_PRIM_SCAL_META_REQUIRE_GENERICS_HPP

#include <stan/math/prim/scal/meta/is_eigen.hpp>
#include <stan/math/prim/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/is_var_or_arithmetic.hpp>
#include <stan/math/prim/scal/meta/is_vector.hpp>
#include <stan/math/prim/scal/meta/is_fvar.hpp>
#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>

#include <type_traits>

namespace stan {

// Alias for most common case
template <bool B>
using bool_constant = std::integral_constant<bool, B>;

/**
 * If condition is true, template is enabled
 */
template <typename Check>
using require_base = std::enable_if_t<Check::value>;

/**
 * If condition is false, template is enabled
 */
template <typename Check>
using require_not = std::enable_if_t<!Check::value>;

/**
 * If all conditions are true, template is enabled
 * Returns a type void if all conditions are true and otherwise fails.
 */
template <class... Checks>
using require_all = std::enable_if_t<math::conjunction<Checks...>::value>;

/**
 * If any condition is true, template is enabled.
 *
 * Returns a type void if any of the conditions are true and otherwise fails.
 */
template <class... Checks>
using require_any = std::enable_if_t<math::disjunction<Checks...>::value>;

/**
 * If all conditions are false, template is enabled.
 *
 * Returns a type void if all of the conditions are false.
 */
template <class... Checks>
using require_all_not = std::enable_if_t<!math::conjunction<Checks...>::value>;

/**
 * If any condition is false, template is enabled.
 *
 * Returns a type void if any of the conditions are false.
 */
template <class... Checks>
using require_any_not = std::enable_if_t<!math::disjunction<Checks...>::value>;

// Enablers for two types of the same value
template <typename T, typename S>
using require_same
    = require_base<std::is_same<std::decay_t<T>, std::decay_t<S>>>;

template <typename T, typename S>
using require_not_same
    = require_not<std::is_same<std::decay_t<T>, std::decay_t<S>>>;

template <typename T, typename... Types>
using require_all_same
    = require_all<std::is_same<std::decay_t<T>, std::decay_t<Types>>...>;

template <typename T, typename... Types>
using not_all_same_type
    = require_all_not<std::is_same<std::decay_t<T>, std::decay_t<Types>>...>;

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
using require_double_or_int = require_base<is_double_or_int<std::decay_t<T>>>;

template <typename T>
using require_not_double_or_int
    = require_not<is_double_or_int<std::decay_t<T>>>;

template <typename... Types>
using require_all_double_or_int
    = require_all<is_double_or_int<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_double_or_int
    = require_any<is_double_or_int<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_double_or_int
    = require_all_not<is_double_or_int<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_double_or_int
    = require_any_not<is_double_or_int<std::decay_t<Types>>...>;

// Checks for arithmetic types
template <typename T>
using require_arithmetic = require_base<std::is_arithmetic<std::decay_t<T>>>;

template <typename T>
using require_not_arithmetic = require_not<std::is_arithmetic<std::decay_t<T>>>;

template <typename... Types>
using require_all_arithmetic
    = require_all<std::is_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_arithmetic
    = require_any<std::is_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_arithmetic
    = require_all_not<std::is_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_arithmetic
    = require_any_not<std::is_arithmetic<std::decay_t<Types>>...>;

// Checks for floating_point types
template <typename T>
using require_floating_point
    = require_base<std::is_floating_point<std::decay_t<T>>>;

template <typename T>
using require_not_floating_point
    = require_not<std::is_floating_point<std::decay_t<T>>>;

template <typename... Types>
using require_all_floating_point
    = require_all<std::is_floating_point<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_floating_point
    = require_any<std::is_floating_point<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_floating_point
    = require_all_not<std::is_floating_point<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_floating_point
    = require_any_not<std::is_floating_point<std::decay_t<Types>>...>;

template <typename T>
using require_var = require_base<is_var<std::decay_t<T>>>;

template <typename T>
using require_not_var = require_not<is_var<std::decay_t<T>>>;

template <typename... Types>
using require_all_var = require_all<is_var<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_var = require_any<is_var<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_var = require_all_not<is_var<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_var = require_any_not<is_var<std::decay_t<Types>>...>;

template <typename T>
using require_var_or_arithmetic
    = require_base<is_var_or_arithmetic<std::decay_t<T>>>;

template <typename T>
using require_not_var_or_arithmetic
    = require_not<is_var_or_arithmetic<std::decay_t<T>>>;

template <typename... Types>
using require_all_var_or_arithmetic
    = require_all<is_var_or_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_var_or_arithmetic
    = require_any<is_var_or_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_var_or_arithmetic
    = require_all_not<is_var_or_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_var_or_arithmetic
    = require_any_not<is_var_or_arithmetic<std::decay_t<Types>>...>;

template <typename T>
using require_fvar = require_base<is_fvar<std::decay_t<T>>>;

template <typename T>
using require_not_fvar = require_not<is_fvar<std::decay_t<T>>>;

template <typename... Types>
using require_all_fvar = require_all<is_fvar<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_fvar = require_any<is_fvar<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_fvar = require_all_not<is_fvar<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_fvar = require_any_not<is_fvar<std::decay_t<Types>>...>;

/**
 * Checks if decayed type is a var or fvar
 * @tparam The type to check
 */
template <typename T>
struct is_var_or_fvar
    : bool_constant<math::disjunction<is_var<std::decay_t<T>>,
                                      is_fvar<std::decay_t<T>>>::value> {};

template <typename T>
using require_var_or_fvar = require_base<is_var_or_fvar<T>>;

template <typename T>
using require_not_var_or_fvar = require_not<is_var_or_fvar<T>>;

template <typename... Types>
using require_all_var_or_fvar = require_all<is_var_or_fvar<Types>...>;

template <typename... Types>
using require_any_var_or_fvar = require_any<is_var_or_fvar<Types>...>;

template <typename... Types>
using require_all_not_var_or_fvar = require_all_not<is_var_or_fvar<Types>...>;

template <typename... Types>
using require_any_not_var_or_fvar = require_any_not<is_var_or_fvar<Types>...>;

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
using require_stan_scalar = require_base<is_stan_scalar<T>>;

template <typename T>
using require_not_stan_scalar = require_not<is_stan_scalar<T>>;

template <typename... Types>
using require_all_stan_scalar = require_all<is_stan_scalar<Types>...>;

template <typename... Types>
using require_any_stan_scalar = require_any<is_stan_scalar<Types>...>;

template <typename... Types>
using require_all_not_stan_scalar = require_all_not<is_stan_scalar<Types>...>;

template <typename... Types>
using require_any_not_stan_scalar = require_any_not<is_stan_scalar<Types>...>;

/**
 * Below are enablers for std_vector based matrices and scalar types of
 * - Double or Int
 * - Arithmetic
 * - Floating Point
 * - Var
 * - Var or Arithmetic
 * - Fvar
 * - Var or Fvar
 */

template <typename T>
using require_std_vector = require_base<is_std_vector<std::decay_t<T>>>;

template <typename T>
using require_not_std_vector = require_not<is_std_vector<std::decay_t<T>>>;

template <typename... Types>
using require_all_std_vector
    = require_all<is_std_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_std_vector
    = require_any<is_std_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_std_vector
    = require_all_not<is_std_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_std_vector
    = require_any_not<is_std_vector<std::decay_t<Types>>...>;

/**
 * Meta struct to check if type is std_vector based and if scalar type passes a
 * conditional test
 * @tparam T The Container to check
 * @tparam CheckType templated type to check value of scalar
 */
template <class Check, template <class> class CheckType>
struct is_std_vector_and_scalar_type
    : bool_constant<
          math::conjunction<is_std_vector<std::decay_t<Check>>,
                            CheckType<scalar_type_decay_t<Check>>>::value> {};

/**
 * Enable if type is std_vector and scalar is double or int
 * @tparam The type to check
 */
template <typename T>
struct is_std_vector_double_or_int
    : is_std_vector_and_scalar_type<T, is_double_or_int> {};

template <typename T>
using require_std_vector_double_or_int
    = require_base<is_std_vector_double_or_int<T>>;

template <typename T>
using require_not_std_vector_double_or_int
    = require_not<is_std_vector_double_or_int<T>>;

template <typename... Types>
using require_all_std_vector_double_or_int
    = require_all<is_std_vector_double_or_int<Types>...>;

template <typename... Types>
using require_any_std_vector_double_or_int
    = require_any<is_std_vector_double_or_int<Types>...>;

template <typename... Types>
using require_all_not_std_vector_double_or_int
    = require_all_not<is_std_vector_double_or_int<Types>...>;

template <typename... Types>
using require_any_not_std_vector_double_or_int
    = require_any_not<is_std_vector_double_or_int<Types>...>;

/**
 * Enable if type is std_vector and scalar is floating point
 * @tparam The type to check
 */
template <typename T>
struct is_std_vector_floating_point
    : is_std_vector_and_scalar_type<T, std::is_floating_point> {};

template <typename T>
using require_std_vector_floating_point
    = require_base<is_std_vector_floating_point<T>>;

template <typename T>
using require_not_std_vector_floating_point
    = require_not<is_std_vector_floating_point<T>>;

template <typename... Types>
using require_all_std_vector_floating_point
    = require_all<is_std_vector_floating_point<Types>...>;

template <typename... Types>
using require_any_std_vector_floating_point
    = require_any<is_std_vector_floating_point<Types>...>;

template <typename... Types>
using require_all_not_std_vector_floating_point
    = require_all_not<is_std_vector_floating_point<Types>...>;

template <typename... Types>
using require_any_not_std_vector_floating_point
    = require_any_not<is_std_vector_floating_point<Types>...>;

/**
 * Enable if type is std_vector and scalar is arithmetic
 * @tparam The type to check
 */
template <typename T>
struct is_std_vector_arithmetic
    : is_std_vector_and_scalar_type<T, std::is_arithmetic> {};

template <typename T>
using require_std_vector_arithmetic = require_base<is_std_vector_arithmetic<T>>;

template <typename T>
using require_not_std_vector_arithmetic
    = require_not<is_std_vector_arithmetic<T>>;

template <typename... Types>
using require_all_std_vector_arithmetic
    = require_all<is_std_vector_arithmetic<Types>...>;

template <typename... Types>
using require_any_std_vector_arithmetic
    = require_any<is_std_vector_arithmetic<Types>...>;

template <typename... Types>
using require_all_not_std_vector_arithmetic
    = require_all_not<is_std_vector_arithmetic<Types>...>;

template <typename... Types>
using require_any_not_std_vector_arithmetic
    = require_any_not<is_std_vector_arithmetic<Types>...>;

/**
 * Enable if type is std_vector and scalar is var
 * @tparam The type to check
 */
template <typename T>
struct is_std_vector_var : is_std_vector_and_scalar_type<T, is_var> {};

template <typename T>
using require_std_vector_var = require_base<is_std_vector_var<T>>;

template <typename T>
using require_not_std_vector_var = require_not<is_std_vector_var<T>>;

template <typename... Types>
using require_all_std_vector_var = require_all<is_std_vector_var<Types>...>;

template <typename... Types>
using require_any_std_vector_var = require_any<is_std_vector_var<Types>...>;

template <typename... Types>
using require_all_not_std_vector_var
    = require_all_not<is_std_vector_var<Types>...>;

template <typename... Types>
using require_any_not_std_vector_var
    = require_any_not<is_std_vector_var<Types>...>;

/**
 * Enable if type is std_vector and scalar is var or arithmetic
 * @tparam The type to check
 */
template <typename T>
struct is_std_vector_var_or_arithmetic
    : is_std_vector_and_scalar_type<T, is_var_or_arithmetic> {};

template <typename T>
using require_std_vector_var_or_arithmetic
    = require_base<is_std_vector_var_or_arithmetic<T>>;

template <typename T>
using require_not_std_vector_var_or_arithmetic
    = require_not<is_std_vector_var_or_arithmetic<T>>;

template <typename... Types>
using require_all_std_vector_var_or_arithmetic
    = require_all<is_std_vector_var_or_arithmetic<Types>...>;

template <typename... Types>
using require_any_std_vector_var_or_arithmetic
    = require_any<is_std_vector_var_or_arithmetic<Types>...>;

template <typename... Types>
using require_all_not_std_vector_var_or_arithmetic
    = require_all_not<is_std_vector_var_or_arithmetic<Types>...>;

template <typename... Types>
using require_any_not_std_vector_var_or_arithmetic
    = require_any_not<is_std_vector_var_or_arithmetic<Types>...>;

/**
 * Enable if type is std_vector and scalar is fvar
 * @tparam The type to check
 */
template <typename T>
struct is_std_vector_fvar : is_std_vector_and_scalar_type<T, is_fvar> {};

template <typename T>
using require_std_vector_fvar = require_base<is_std_vector_fvar<T>>;

template <typename T>
using require_not_std_vector_fvar = require_not<is_std_vector_fvar<T>>;

template <typename... Types>
using require_all_std_vector_fvar = require_all<is_std_vector_fvar<Types>...>;

template <typename... Types>
using require_any_std_vector_fvar = require_any<is_std_vector_fvar<Types>...>;

template <typename... Types>
using require_all_not_std_vector_fvar
    = require_all_not<is_std_vector_fvar<Types>...>;

template <typename... Types>
using require_any_not_std_vector_fvar
    = require_any_not<is_std_vector_fvar<Types>...>;

/**
 * Enable if type is std_vector and scalar is var or fvar
 * @tparam The type to check
 */
template <typename T>
struct is_std_vector_var_or_fvar
    : is_std_vector_and_scalar_type<T, is_var_or_fvar> {};

template <typename T>
using require_std_vector_var_or_fvar
    = require_base<is_std_vector_var_or_fvar<T>>;

template <typename T>
using require_not_std_vector_var_or_fvar
    = require_not<is_std_vector_var_or_fvar<T>>;

template <typename... Types>
using require_all_std_vector_var_or_fvar
    = require_all<is_std_vector_var_or_fvar<Types>...>;

template <typename... Types>
using require_any_std_vector_var_or_fvar
    = require_any<is_std_vector_var_or_fvar<Types>...>;

template <typename... Types>
using require_all_not_std_vector_var_or_fvar
    = require_all_not<is_std_vector_var_or_fvar<Types>...>;

template <typename... Types>
using require_any_not_std_vector_var_or_fvar
    = require_any_not<is_std_vector_var_or_fvar<Types>...>;

/**
 * Below are enablers for vector based matrices and scalar types of
 * - Double or Int
 * - Arithmetic
 * - Floating Point
 * - Var
 * - Var or Arithmetic
 * - Fvar
 * - Var or Fvar
 */

template <typename T>
using require_vector = require_base<is_vector<std::decay_t<T>>>;

template <typename T>
using require_not_vector = require_not<is_vector<std::decay_t<T>>>;

template <typename... Types>
using require_all_vector = require_all<is_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_vector = require_any<is_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_vector
    = require_all_not<is_vector<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_vector
    = require_any_not<is_vector<std::decay_t<Types>>...>;

/**
 * Meta struct to check if type is vector based and if scalar type passes a
 * conditional test
 * @tparam T The Container to check
 * @tparam CheckType templated type to check value of scalar
 */
template <class Check, template <class> class CheckType>
struct is_vector_and_scalar_type
    : bool_constant<
          math::conjunction<is_vector<std::decay_t<Check>>,
                            CheckType<scalar_type_decay_t<Check>>>::value> {};

/**
 * Enable if type is vector and scalar is double or int
 * @tparam The type to check
 */
template <typename T>
struct is_vector_double_or_int
    : is_vector_and_scalar_type<T, is_double_or_int> {};

template <typename T>
using require_vector_double_or_int = require_base<is_vector_double_or_int<T>>;

template <typename T>
using require_not_vector_double_or_int
    = require_not<is_vector_double_or_int<T>>;

template <typename... Types>
using require_all_vector_double_or_int
    = require_all<is_vector_double_or_int<Types>...>;

template <typename... Types>
using require_any_vector_double_or_int
    = require_any<is_vector_double_or_int<Types>...>;

template <typename... Types>
using require_all_not_vector_double_or_int
    = require_all_not<is_vector_double_or_int<Types>...>;

template <typename... Types>
using require_any_not_vector_double_or_int
    = require_any_not<is_vector_double_or_int<Types>...>;

/**
 * Enable if type is vector and scalar is floating point
 * @tparam The type to check
 */
template <typename T>
struct is_vector_floating_point
    : is_vector_and_scalar_type<T, std::is_floating_point> {};

template <typename T>
using require_vector_floating_point = require_base<is_vector_floating_point<T>>;

template <typename T>
using require_not_vector_floating_point
    = require_not<is_vector_floating_point<T>>;

template <typename... Types>
using require_all_vector_floating_point
    = require_all<is_vector_floating_point<Types>...>;

template <typename... Types>
using require_any_vector_floating_point
    = require_any<is_vector_floating_point<Types>...>;

template <typename... Types>
using require_all_not_vector_floating_point
    = require_all_not<is_vector_floating_point<Types>...>;

template <typename... Types>
using require_any_not_vector_floating_point
    = require_any_not<is_vector_floating_point<Types>...>;

/**
 * Enable if type is vector and scalar is arithmetic
 * @tparam The type to check
 */
template <typename T>
struct is_vector_arithmetic : is_vector_and_scalar_type<T, std::is_arithmetic> {
};

template <typename T>
using require_vector_arithmetic = require_base<is_vector_arithmetic<T>>;

template <typename T>
using require_not_vector_arithmetic = require_not<is_vector_arithmetic<T>>;

template <typename... Types>
using require_all_vector_arithmetic
    = require_all<is_vector_arithmetic<Types>...>;

template <typename... Types>
using require_any_vector_arithmetic
    = require_any<is_vector_arithmetic<Types>...>;

template <typename... Types>
using require_all_not_vector_arithmetic
    = require_all_not<is_vector_arithmetic<Types>...>;

template <typename... Types>
using require_any_not_vector_arithmetic
    = require_any_not<is_vector_arithmetic<Types>...>;

/**
 * Enable if type is vector and scalar is var
 * @tparam The type to check
 */
template <typename T>
struct is_vector_var : is_vector_and_scalar_type<T, is_var> {};

template <typename T>
using require_vector_var = require_base<is_vector_var<T>>;

template <typename T>
using require_not_vector_var = require_not<is_vector_var<T>>;

template <typename... Types>
using require_all_vector_var = require_all<is_vector_var<Types>...>;

template <typename... Types>
using require_any_vector_var = require_any<is_vector_var<Types>...>;

template <typename... Types>
using require_all_not_vector_var = require_all_not<is_vector_var<Types>...>;

template <typename... Types>
using require_any_not_vector_var = require_any_not<is_vector_var<Types>...>;

/**
 * Enable if type is vector and scalar is var or arithmetic
 * @tparam The type to check
 */
template <typename T>
struct is_vector_var_or_arithmetic
    : is_vector_and_scalar_type<T, is_var_or_arithmetic> {};

template <typename T>
using require_vector_var_or_arithmetic
    = require_base<is_vector_var_or_arithmetic<T>>;

template <typename T>
using require_not_vector_var_or_arithmetic
    = require_not<is_vector_var_or_arithmetic<T>>;

template <typename... Types>
using require_all_vector_var_or_arithmetic
    = require_all<is_vector_var_or_arithmetic<Types>...>;

template <typename... Types>
using require_any_vector_var_or_arithmetic
    = require_any<is_vector_var_or_arithmetic<Types>...>;

template <typename... Types>
using require_all_not_vector_var_or_arithmetic
    = require_all_not<is_vector_var_or_arithmetic<Types>...>;

template <typename... Types>
using require_any_not_vector_var_or_arithmetic
    = require_any_not<is_vector_var_or_arithmetic<Types>...>;

/**
 * Enable if type is vector and scalar is fvar
 * @tparam The type to check
 */
template <typename T>
struct is_vector_fvar : is_vector_and_scalar_type<T, is_fvar> {};

template <typename T>
using require_vector_fvar = require_base<is_vector_fvar<T>>;

template <typename T>
using require_not_vector_fvar = require_not<is_vector_fvar<T>>;

template <typename... Types>
using require_all_vector_fvar = require_all<is_vector_fvar<Types>...>;

template <typename... Types>
using require_any_vector_fvar = require_any<is_vector_fvar<Types>...>;

template <typename... Types>
using require_all_not_vector_fvar = require_all_not<is_vector_fvar<Types>...>;

template <typename... Types>
using require_any_not_vector_fvar = require_any_not<is_vector_fvar<Types>...>;

/**
 * Enable if type is vector and scalar is var or fvar
 * @tparam The type to check
 */
template <typename T>
struct is_vector_var_or_fvar : is_vector_and_scalar_type<T, is_var_or_fvar> {};

template <typename T>
using require_vector_var_or_fvar = require_base<is_vector_var_or_fvar<T>>;

template <typename T>
using require_not_vector_var_or_fvar = require_not<is_vector_var_or_fvar<T>>;

template <typename... Types>
using require_all_vector_var_or_fvar
    = require_all<is_vector_var_or_fvar<Types>...>;

template <typename... Types>
using require_any_vector_var_or_fvar
    = require_any<is_vector_var_or_fvar<Types>...>;

template <typename... Types>
using require_all_not_vector_var_or_fvar
    = require_all_not<is_vector_var_or_fvar<Types>...>;

template <typename... Types>
using require_any_not_vector_var_or_fvar
    = require_any_not<is_vector_var_or_fvar<Types>...>;

/**
 * Below are enablers for eigen based matrices and scalar types of
 * - Double or Int
 * - Arithmetic
 * - Floating Point
 * - Var
 * - Var or Arithmetic
 * - Fvar
 * - Var or Fvar
 */

template <typename T>
using require_eigen = require_base<is_eigen<std::decay_t<T>>>;

template <typename T>
using require_not_eigen = require_not<is_eigen<std::decay_t<T>>>;

template <typename... Types>
using require_all_eigen = require_all<is_eigen<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_eigen = require_any<is_eigen<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_eigen = require_all_not<is_eigen<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_eigen = require_any_not<is_eigen<std::decay_t<Types>>...>;

/**
 * Meta struct to check if type is eigen based and if scalar type passes a
 * conditional test
 * @tparam T The Container to check
 * @tparam CheckType templated type to check value of scalar
 */
template <class Check, template <class> class CheckType>
struct is_eigen_and_scalar_type
    : bool_constant<
          math::conjunction<is_eigen<std::decay_t<Check>>,
                            CheckType<scalar_type_decay_t<Check>>>::value> {};

/**
 * Enable if type is eigen and scalar is double or int
 * @tparam The type to check
 */
template <typename T>
struct is_eigen_double_or_int : is_eigen_and_scalar_type<T, is_double_or_int> {
};

template <typename T>
using require_eigen_double_or_int = require_base<is_eigen_double_or_int<T>>;

template <typename T>
using require_not_eigen_double_or_int = require_not<is_eigen_double_or_int<T>>;

template <typename... Types>
using require_all_eigen_double_or_int
    = require_all<is_eigen_double_or_int<Types>...>;

template <typename... Types>
using require_any_eigen_double_or_int
    = require_any<is_eigen_double_or_int<Types>...>;

template <typename... Types>
using require_all_not_eigen_double_or_int
    = require_all_not<is_eigen_double_or_int<Types>...>;

template <typename... Types>
using require_any_not_eigen_double_or_int
    = require_any_not<is_eigen_double_or_int<Types>...>;

/**
 * Enable if type is eigen and scalar is floating point
 * @tparam The type to check
 */
template <typename T>
struct is_eigen_floating_point
    : is_eigen_and_scalar_type<T, std::is_floating_point> {};

template <typename T>
using require_eigen_floating_point = require_base<is_eigen_floating_point<T>>;

template <typename T>
using require_not_eigen_floating_point
    = require_not<is_eigen_floating_point<T>>;

template <typename... Types>
using require_all_eigen_floating_point
    = require_all<is_eigen_floating_point<Types>...>;

template <typename... Types>
using require_any_eigen_floating_point
    = require_any<is_eigen_floating_point<Types>...>;

template <typename... Types>
using require_all_not_eigen_floating_point
    = require_all_not<is_eigen_floating_point<Types>...>;

template <typename... Types>
using require_any_not_eigen_floating_point
    = require_any_not<is_eigen_floating_point<Types>...>;

/**
 * Enable if type is eigen and scalar is arithmetic
 * @tparam The type to check
 */
template <typename T>
struct is_eigen_arithmetic : is_eigen_and_scalar_type<T, std::is_arithmetic> {};

template <typename T>
using require_eigen_arithmetic = require_base<is_eigen_arithmetic<T>>;

template <typename T>
using require_not_eigen_arithmetic = require_not<is_eigen_arithmetic<T>>;

template <typename... Types>
using require_all_eigen_arithmetic = require_all<is_eigen_arithmetic<Types>...>;

template <typename... Types>
using require_any_eigen_arithmetic = require_any<is_eigen_arithmetic<Types>...>;

template <typename... Types>
using require_all_not_eigen_arithmetic
    = require_all_not<is_eigen_arithmetic<Types>...>;

template <typename... Types>
using require_any_not_eigen_arithmetic
    = require_any_not<is_eigen_arithmetic<Types>...>;

/**
 * Enable if type is eigen and scalar is var
 * @tparam The type to check
 */
template <typename T>
struct is_eigen_var : is_eigen_and_scalar_type<T, is_var> {};

template <typename T>
using require_eigen_var = require_base<is_eigen_var<T>>;

template <typename T>
using require_not_eigen_var = require_not<is_eigen_var<T>>;

template <typename... Types>
using require_all_eigen_var = require_all<is_eigen_var<Types>...>;

template <typename... Types>
using require_any_eigen_var = require_any<is_eigen_var<Types>...>;

template <typename... Types>
using require_all_not_eigen_var = require_all_not<is_eigen_var<Types>...>;

template <typename... Types>
using require_any_not_eigen_var = require_any_not<is_eigen_var<Types>...>;

/**
 * Enable if type is eigen and scalar is var or arithmetic
 * @tparam The type to check
 */
template <typename T>
struct is_eigen_var_or_arithmetic
    : is_eigen_and_scalar_type<T, is_var_or_arithmetic> {};

template <typename T>
using require_eigen_var_or_arithmetic
    = require_base<is_eigen_var_or_arithmetic<T>>;

template <typename T>
using require_not_eigen_var_or_arithmetic
    = require_not<is_eigen_var_or_arithmetic<T>>;

template <typename... Types>
using require_all_eigen_var_or_arithmetic
    = require_all<is_eigen_var_or_arithmetic<Types>...>;

template <typename... Types>
using require_any_eigen_var_or_arithmetic
    = require_any<is_eigen_var_or_arithmetic<Types>...>;

template <typename... Types>
using require_all_not_eigen_var_or_arithmetic
    = require_all_not<is_eigen_var_or_arithmetic<Types>...>;

template <typename... Types>
using require_any_not_eigen_var_or_arithmetic
    = require_any_not<is_eigen_var_or_arithmetic<Types>...>;

/**
 * Enable if type is eigen and scalar is fvar
 * @tparam The type to check
 */
template <typename T>
struct is_eigen_fvar : is_eigen_and_scalar_type<T, is_fvar> {};

template <typename T>
using require_eigen_fvar = require_base<is_eigen_fvar<T>>;

template <typename T>
using require_not_eigen_fvar = require_not<is_eigen_fvar<T>>;

template <typename... Types>
using require_all_eigen_fvar = require_all<is_eigen_fvar<Types>...>;

template <typename... Types>
using require_any_eigen_fvar = require_any<is_eigen_fvar<Types>...>;

template <typename... Types>
using require_all_not_eigen_fvar = require_all_not<is_eigen_fvar<Types>...>;

template <typename... Types>
using require_any_not_eigen_fvar = require_any_not<is_eigen_fvar<Types>...>;

/**
 * Enable if type is eigen and scalar is var or fvar
 * @tparam The type to check
 */
template <typename T>
struct is_eigen_var_or_fvar : is_eigen_and_scalar_type<T, is_var_or_fvar> {};

template <typename T>
using require_eigen_var_or_fvar = require_base<is_eigen_var_or_fvar<T>>;

template <typename T>
using require_not_eigen_var_or_fvar = require_not<is_eigen_var_or_fvar<T>>;

template <typename... Types>
using require_all_eigen_var_or_fvar
    = require_all<is_eigen_var_or_fvar<Types>...>;

template <typename... Types>
using require_any_eigen_var_or_fvar
    = require_any<is_eigen_var_or_fvar<Types>...>;

template <typename... Types>
using require_all_not_eigen_var_or_fvar
    = require_all_not<is_eigen_var_or_fvar<Types>...>;

template <typename... Types>
using require_any_not_eigen_var_or_fvar
    = require_any_not<is_eigen_var_or_fvar<Types>...>;

}  // namespace stan
#endif
