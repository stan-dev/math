#ifndef STAN_MATH_PRIM_SCAL_META_REQUIRE_GENERICS_HPP
#define STAN_MATH_PRIM_SCAL_META_REQUIRE_GENERICS_HPP

#include <stan/math/prim/scal/meta/is_eigen.hpp>
#include <stan/math/prim/scal/meta/is_vector.hpp>
#include <stan/math/prim/scal/meta/is_var.hpp>
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
using require_double_or_int = require_base<is_double_or_int<T>>;

template <typename T>
using require_not_double_or_int = require_not<is_double_or_int<T>>;

template <typename... Types>
using require_all_double_or_int = require_all<is_double_or_int<Types>...>;

template <typename... Types>
using require_any_double_or_int = require_any<is_double_or_int<Types>...>;

template <typename... Types>
using require_all_not_double_or_int
    = require_all_not<is_double_or_int<Types>...>;

template <typename... Types>
using require_any_not_double_or_int
    = require_any_not<is_double_or_int<Types>...>;

/**
 * Meta struct to check if type is eigen if scalar type passes a conditional
 * test
 * @tparam The type to check
 */
template <class Check, template <class> class CheckType>
struct is_eigen_and_scalar_type
    : bool_constant<
          math::conjunction<is_eigen<std::decay_t<Check>>,
                            CheckType<scalar_type_decay_t<Check>>>::value> {};

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
using require_not_all_eigen_arithmetic
    = require_all_not<is_eigen_arithmetic<Types>...>;

template <typename... Types>
using require_not_any_eigen_arithmetic
    = require_any_not<is_eigen_arithmetic<Types>...>;

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
using require_not_all_eigen_floating_point
    = require_all_not<is_eigen_floating_point<Types>...>;

template <typename... Types>
using require_not_any_eigen_floating_point
    = require_any_not<is_eigen_floating_point<Types>...>;

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
using require_not_all_eigen_var = require_all_not<is_eigen_var<Types>...>;

template <typename... Types>
using require_not_any_eigen_var = require_any_not<is_eigen_var<Types>...>;

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
using require_not_all_eigen_fvar = require_all_not<is_eigen_fvar<Types>...>;

template <typename... Types>
using require_not_any_eigen_fvar = require_any_not<is_eigen_fvar<Types>...>;

template <class Check, template <class> class CheckType>
struct is_std_vector_and_scalar_type
    : bool_constant<math::conjunction<
          is_std_vector<std::decay_t<Check>>,
          CheckType<scalar_type_decay_t<std::decay_t<Check>>>>::value> {};

// Check if a type is std_vector with arithmetic underlying type
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
using require_not_all_std_vector_arithmetic
    = require_all_not<is_std_vector_arithmetic<Types>...>;

template <typename... Types>
using require_not_any_std_vector_arithmetic
    = require_any_not<is_std_vector_arithmetic<Types>...>;

// Check if a type is std_vector with floating point underlying type
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
using require_not_all_std_vector_floating_point
    = require_all_not<is_std_vector_floating_point<Types>...>;

template <typename... Types>
using require_not_any_std_vector_floating_point
    = require_any_not<is_std_vector_floating_point<Types>...>;

// Check if a type is std_vector with floatingpoint underlying type
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
using require_not_all_std_vector_var
    = require_all_not<is_std_vector_var<Types>...>;

template <typename... Types>
using require_not_any_std_vector_var
    = require_any_not<is_std_vector_var<Types>...>;

// Check if a type is std_vector with floatingpoint underlying type
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
using require_not_all_std_vector_fvar
    = require_all_not<is_std_vector_fvar<Types>...>;

template <typename... Types>
using require_not_any_std_vector_fvar
    = require_any_not<is_std_vector_fvar<Types>...>;

}  // namespace stan
#endif
