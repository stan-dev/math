#ifndef STAN_MATH_PRIM_SCAL_META_REQUIRE_GENERICS_HPP
#define STAN_MATH_PRIM_SCAL_META_REQUIRE_GENERICS_HPP

#include <stan/math/prim/scal/meta/bool_constant.hpp>
#include <stan/math/prim/scal/meta/is_eigen.hpp>
#include <stan/math/prim/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/is_var_or_arithmetic.hpp>
#include <stan/math/prim/scal/meta/is_vector.hpp>
#include <stan/math/prim/scal/meta/is_fvar.hpp>
#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <stan/math/prim/scal/meta/value_type.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>

#include <type_traits>

namespace stan {

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
using require_not_all_t
    = std::enable_if_t<!math::conjunction<Checks...>::value>;

/**
 * If any condition is false, template is enabled.
 *
 * Returns a type void if any of the conditions are false.
 */
template <class... Checks>
using require_not_any_t
    = std::enable_if_t<!math::disjunction<Checks...>::value>;

// Enablers for two types of the same value
template <typename T, typename S>
using require_same = require_t<std::is_same<std::decay_t<T>, std::decay_t<S>>>;

template <typename T, typename S>
using require_not_same
    = require_not_t<std::is_same<std::decay_t<T>, std::decay_t<S>>>;

template <typename T, typename... Types>
using require_all_same
    = require_all_t<std::is_same<std::decay_t<T>, std::decay_t<Types>>...>;

template <typename T, typename... Types>
using require_not_all_same
    = require_not_all_t<std::is_same<std::decay_t<T>, std::decay_t<Types>>...>;

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
using require_double_or_int = require_t<is_double_or_int<std::decay_t<T>>>;

template <typename T>
using require_not_double_or_int
    = require_not_t<is_double_or_int<std::decay_t<T>>>;

template <typename... Types>
using require_all_double_or_int
    = require_all_t<is_double_or_int<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_double_or_int
    = require_any_t<is_double_or_int<std::decay_t<Types>>...>;

template <typename... Types>
using require_not_all_double_or_int
    = require_not_all_t<is_double_or_int<std::decay_t<Types>>...>;

template <typename... Types>
using require_not_any_double_or_int
    = require_not_any_t<is_double_or_int<std::decay_t<Types>>...>;

// Checks for arithmetic types
template <typename T>
using require_arithmetic = require_t<std::is_arithmetic<std::decay_t<T>>>;

template <typename T>
using require_not_arithmetic
    = require_not_t<std::is_arithmetic<std::decay_t<T>>>;

template <typename... Types>
using require_all_arithmetic
    = require_all_t<std::is_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_arithmetic
    = require_any_t<std::is_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_not_all_arithmetic
    = require_not_all_t<std::is_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_not_any_arithmetic
    = require_not_any_t<std::is_arithmetic<std::decay_t<Types>>...>;

// Checks for floating_point types
template <typename T>
using require_floating_point
    = require_t<std::is_floating_point<std::decay_t<T>>>;

template <typename T>
using require_not_floating_point
    = require_not_t<std::is_floating_point<std::decay_t<T>>>;

template <typename... Types>
using require_all_floating_point
    = require_all_t<std::is_floating_point<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_floating_point
    = require_any_t<std::is_floating_point<std::decay_t<Types>>...>;

template <typename... Types>
using require_not_all_floating_point
    = require_not_all_t<std::is_floating_point<std::decay_t<Types>>...>;

template <typename... Types>
using require_not_any_floating_point
    = require_not_any_t<std::is_floating_point<std::decay_t<Types>>...>;

template <typename T>
using require_var = require_t<is_var<std::decay_t<T>>>;

template <typename T>
using require_not_var = require_not_t<is_var<std::decay_t<T>>>;

template <typename... Types>
using require_all_var = require_all_t<is_var<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_var = require_any_t<is_var<std::decay_t<Types>>...>;

template <typename... Types>
using require_not_all_var = require_not_all_t<is_var<std::decay_t<Types>>...>;

template <typename... Types>
using require_not_any_var = require_not_any_t<is_var<std::decay_t<Types>>...>;

template <typename T>
using require_var_or_arithmetic
    = require_t<is_var_or_arithmetic<std::decay_t<T>>>;

template <typename T>
using require_not_var_or_arithmetic
    = require_not_t<is_var_or_arithmetic<std::decay_t<T>>>;

template <typename... Types>
using require_all_var_or_arithmetic
    = require_all_t<is_var_or_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_var_or_arithmetic
    = require_any_t<is_var_or_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_not_all_var_or_arithmetic
    = require_not_all_t<is_var_or_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_not_any_var_or_arithmetic
    = require_not_any_t<is_var_or_arithmetic<std::decay_t<Types>>...>;

template <typename T>
using require_fvar = require_t<is_fvar<std::decay_t<T>>>;

template <typename T>
using require_not_fvar = require_not_t<is_fvar<std::decay_t<T>>>;

template <typename... Types>
using require_all_fvar = require_all_t<is_fvar<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_fvar = require_any_t<is_fvar<std::decay_t<Types>>...>;

template <typename... Types>
using require_not_all_fvar = require_not_all_t<is_fvar<std::decay_t<Types>>...>;

template <typename... Types>
using require_not_any_fvar = require_not_any_t<is_fvar<std::decay_t<Types>>...>;

/**
 * Checks if decayed type is a var or fvar
 * @tparam The type to check
 */
template <typename T>
struct is_var_or_fvar
    : bool_constant<math::disjunction<is_var<std::decay_t<T>>,
                                      is_fvar<std::decay_t<T>>>::value> {};

template <typename T>
using require_var_or_fvar = require_t<is_var_or_fvar<T>>;

template <typename T>
using require_not_var_or_fvar = require_not_t<is_var_or_fvar<T>>;

template <typename... Types>
using require_all_var_or_fvar = require_all_t<is_var_or_fvar<Types>...>;

template <typename... Types>
using require_any_var_or_fvar = require_any_t<is_var_or_fvar<Types>...>;

template <typename... Types>
using require_not_all_var_or_fvar = require_not_all_t<is_var_or_fvar<Types>...>;

template <typename... Types>
using require_not_any_var_or_fvar = require_not_any_t<is_var_or_fvar<Types>...>;

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
using require_stan_scalar = require_t<is_stan_scalar<T>>;

template <typename T>
using require_not_stan_scalar = require_not_t<is_stan_scalar<T>>;

template <typename... Types>
using require_all_stan_scalar = require_all_t<is_stan_scalar<Types>...>;

template <typename... Types>
using require_any_stan_scalar = require_any_t<is_stan_scalar<Types>...>;

template <typename... Types>
using require_not_all_stan_scalar = require_not_all_t<is_stan_scalar<Types>...>;

template <typename... Types>
using require_not_any_stan_scalar = require_not_any_t<is_stan_scalar<Types>...>;

/**
 * Requires for containers
 */

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
using container_type_check_base = bool_constant<
    math::conjunction<ContainerCheck<std::decay_t<Check>>...,
                      TypeCheck<value_type_t<Check>>...>::value>;

/**
 * Check a templated type to see if it and it's inner type pass a condiational
 * test.
 * @tparam ContainerCheck Templated struct or alias that wraps a static constant
 * value called type. Used to check the container satisfies a particular type
 * check. used like template <typename T, require_container_t<is_std_vector,
 * is_var, T>...>
 */
template <template <class...> class ContainerCheck,
          template <class...> class TypeCheck, class... Check>
using require_container_t
    = require_t<container_type_check_base<ContainerCheck, TypeCheck, Check...>>;

/**
 * Std vector
 */
template <template <class...> class TypeCheck, class... Check>
struct is_std_vector_check
    : container_type_check_base<is_std_vector, TypeCheck, Check...> {};

template <template <class...> class TypeCheck, class... Check>
using require_std_vector_t
    = require_t<is_std_vector_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_std_vector_t
    = require_not_t<is_std_vector_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_std_vector_t
    = require_any_t<is_std_vector_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_not_any_std_vector_t
    = require_not_any_t<is_std_vector_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_std_vector_t
    = require_all_t<is_std_vector_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_not_all_std_vector_t
    = require_not_all_t<is_std_vector_check<TypeCheck, Check>...>;

/**
 * Std vector
 */
template <template <class...> class TypeCheck, class... Check>
struct is_vector_check
    : container_type_check_base<is_vector, TypeCheck, Check...> {};

template <template <class...> class TypeCheck, class... Check>
using require_vector_t = require_t<is_vector_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_vector_t
    = require_not_t<is_vector_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_vector_t
    = require_any_t<is_vector_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_not_any_vector_t
    = require_not_any_t<is_vector_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_vector_t
    = require_all_t<is_vector_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_not_all_vector_t
    = require_not_all_t<is_vector_check<TypeCheck, Check>...>;

/**
 * Std vector
 */
template <template <class...> class TypeCheck, class... Check>
struct is_vector_like_check
    : container_type_check_base<is_vector_like, TypeCheck, Check...> {};

template <template <class...> class TypeCheck, class... Check>
using require_vector_like_t
    = require_t<is_vector_like_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_vector_like_t
    = require_not_t<is_vector_like_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_vector_like_t
    = require_any_t<is_vector_like_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_not_any_vector_like_t
    = require_not_any_t<is_vector_like_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_vector_like_t
    = require_all_t<is_vector_like_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_not_all_vector_like_t
    = require_not_all_t<is_vector_like_check<TypeCheck, Check>...>;

/**
 * Std vector
 */
template <template <class...> class TypeCheck, class... Check>
struct is_eigen_check
    : container_type_check_base<is_eigen, TypeCheck, Check...> {};

template <template <class...> class TypeCheck, class... Check>
using require_eigen_t = require_t<is_eigen_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_eigen_t = require_not_t<is_eigen_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_eigen_t = require_any_t<is_eigen_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_not_any_eigen_t
    = require_not_any_t<is_eigen_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_eigen_t = require_all_t<is_eigen_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_not_all_eigen_t
    = require_not_all_t<is_eigen_check<TypeCheck, Check>...>;

/**
 * Std vector
 */
template <template <class...> class TypeCheck, class... Check>
struct is_eigen_vector_check
    : container_type_check_base<is_eigen_vector, TypeCheck, Check...> {};

template <template <class...> class TypeCheck, class... Check>
using require_eigen_vector_t
    = require_t<is_eigen_vector_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_eigen_vector_t
    = require_not_t<is_eigen_vector_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_eigen_vector_t
    = require_any_t<is_eigen_vector_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_not_any_eigen_vector_t
    = require_not_any_t<is_eigen_vector_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_eigen_vector_t
    = require_all_t<is_eigen_vector_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_not_all_eigen_vector_t
    = require_not_all_t<is_eigen_vector_check<TypeCheck, Check>...>;

}  // namespace stan
#endif
