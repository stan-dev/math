#ifndef STAN_MATH_PRIM_SCAL_META_ENABLE_IF_GENERICS_HPP
#define STAN_MATH_PRIM_SCAL_META_ENABLE_IF_GENERICS_HPP

#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>
#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <stan/math/prim/scal/meta/is_var.hpp>
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
using enable_if_not = typename std::enable_if_t<!Check::value>;

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
using enable_if_all_not = std::enable_if_t<!math::conjunction<Checks...>::value>;

/**
 * Template metaprogram for deducing if a method should be remove a
 * a function from overload resolution.
 *
 * Returns a type void if any of the conditions are false.
 */
template <class... Checks>
using enable_if_any_not = std::enable_if_t<!math::disjunction<Checks...>::value>;

// Check whether the decayed types are the same
template <typename T, typename S>
struct is_same_decay
    : std::integral_constant<bool,
                             std::is_same<std::decay_t<T>, std::decay_t<S>>::value> {};

template <typename T, typename S>
using enable_if_same
    = std::enable_if_t<is_same_decay<T, S>::value>;

template <typename T, typename S>
using enable_if_not_same = enable_if_not<is_same_decay<T, S>>;

template <typename T, typename... Types>
using enable_if_all_same = enable_if_all<is_same_decay<T, Types>...>;

template <typename T, typename... Types>
using enable_if_all_not_same = enable_if_all_not<is_same_decay<T, Types>...>;



// Checks decayed (non-const/ref'd) type is arithmetic
template <typename T>
struct is_arithmetic_decay
    : std::integral_constant<bool,
                             std::is_arithmetic<std::decay_t<T>>::value> {};


template <typename T>
using enable_if_arithmetic = std::enable_if_t<is_arithmetic_decay<T>::value>;

template <typename T>
using enable_if_not_arithmetic = enable_if_not<is_arithmetic_decay<T>>;

template <typename... Types>
using enable_if_all_arithmetic = enable_if_all<is_arithmetic_decay<Types>...>;

template <typename... Types>
using enable_if_any_arithmetic = enable_if_any<is_arithmetic_decay<Types>...>;

template <typename... Types>
using enable_if_all_not_arithmetic = enable_if_all_not<is_arithmetic_decay<Types>...>;

template <typename... Types>
using enable_if_any_not_arithmetic = enable_if_any_not<is_arithmetic_decay<Types>...>;

// Check if a type contains or is arithmetic
template <typename T>
struct is_contains_arithmetic
    : std::integral_constant<bool,std::is_arithmetic<scalar_type_t<std::decay_t<T>>>::value> {};

template <typename T>
using enable_if_contains_arithmetic = std::enable_if_t<is_contains_arithmetic<T>::value>;

template <typename T>
using enable_if_not_contains_arithmetic = enable_if_not<is_contains_arithmetic<T>>;

template <typename... Types>
using enable_if_all_contains_arithmetic = enable_if_all<is_contains_arithmetic<Types>...>;

template <typename... Types>
using enable_if_any_contains_arithmetic = enable_if_any<is_contains_arithmetic<Types>...>;

template <typename... Types>
using enable_if_all_not_contains_arithmetic = enable_if_all_not<is_contains_arithmetic<Types>...>;

template <typename... Types>
using enable_if_any_not_contains_arithmetic = enable_if_any_not<is_contains_arithmetic<Types>...>;


// Checks whether the type is floating_point
template <typename T>
struct is_fp_decay
    : std::integral_constant<bool,
                             std::is_floating_point<std::decay_t<T>>::value> {};

template <typename T>
using enable_if_floating_point = std::enable_if_t<is_fp_decay<T>::value>;

template <typename T>
using enable_if_not_floating_point = enable_if_not<is_fp_decay<T>>;

template <typename... Types>
using enable_if_all_floating_point = enable_if_all<is_fp_decay<Types>...>;

template <typename... Types>
using enable_if_any_floating_point = enable_if_any<is_fp_decay<Types>...>;

template <typename... Types>
using enable_if_all_not_floating_point = enable_if_all_not<is_fp_decay<Types>...>;

template <typename... Types>
using enable_if_any_not_floating_point = enable_if_any_not<is_fp_decay<Types>...>;

// Check if a type contains or is floating point
template <typename T>
struct is_contains_floating_point
    : std::integral_constant<bool,std::is_floating_point<scalar_type_t<std::decay_t<T>>>::value> {};

template <typename T>
using enable_if_contains_floating_point = std::enable_if_t<is_contains_floating_point<T>::value>;

template <typename T>
using enable_if_not_contains_floating_point = enable_if_not<is_contains_floating_point<T>>;

template <typename... Types>
using enable_if_all_contains_floating_point = enable_if_all<is_contains_floating_point<Types>...>;

template <typename... Types>
using enable_if_any_contains_floating_point = enable_if_any<is_contains_floating_point<Types>...>;

template <typename... Types>
using enable_if_all_not_contains_floating_point = enable_if_all_not<is_contains_floating_point<Types>...>;

template <typename... Types>
using enable_if_any_not_contains_floating_point = enable_if_any_not<is_contains_floating_point<Types>...>;

// Check if type is a var

template <typename T>
using enable_if_var = std::enable_if_t<is_var<std::decay<T>>::value>;

template <typename T>
using enable_if_not_var = enable_if_not<is_var<std::decay<T>>>;

template <typename... Types>
using enable_if_all_var = enable_if_all<is_var<std::decay<Types>>...>;

template <typename... Types>
using enable_if_any_var = enable_if_any<is_var<std::decay<Types>>...>;

template <typename... Types>
using enable_if_all_not_var = enable_if_all_not<is_var<std::decay<Types>>...>;

template <typename... Types>
using enable_if_any_not_var = enable_if_any_not<is_var<std::decay<Types>>...>;


// Check if type contains or is a var
template <typename T>
struct is_contains_var
    : std::integral_constant<bool,
                             is_var<scalar_type_t<std::decay_t<T>>>::value> {};

template <typename T>
using enable_if_contains_var = std::enable_if_t<is_contains_var<T>::value>;

template <typename T>
using enable_if_not_contains_var = enable_if_not<is_contains_var<T>>;

template <typename... Types>
using enable_if_all_contains_var = enable_if_all<is_contains_var<Types>...>;

template <typename... Types>
using enable_if_any_contains_var = enable_if_any<is_contains_var<Types>...>;

template <typename... Types>
using enable_if_all_not_contains_var = enable_if_all_not<is_contains_var<Types>...>;

template <typename... Types>
using enable_if_any_not_contains_var = enable_if_any_not<is_contains_var<Types>...>;

// Check if type is a fvar

template <typename T>
using enable_if_fvar = std::enable_if_t<is_fvar<std::decay<T>>::value>;

template <typename T>
using enable_if_not_fvar = enable_if_not<is_fvar<std::decay<T>>>;

template <typename... Types>
using enable_if_all_fvar = enable_if_all<is_fvar<std::decay<Types>>...>;

template <typename... Types>
using enable_if_any_fvar = enable_if_any<is_fvar<std::decay<Types>>...>;

template <typename... Types>
using enable_if_all_not_fvar = enable_if_all_not<is_fvar<std::decay<Types>>...>;

template <typename... Types>
using enable_if_any_not_fvar = enable_if_any_not<is_fvar<std::decay<Types>>...>;


// Check if type contains or is a fvar
template <typename T>
struct is_contains_fvar
    : std::integral_constant<bool,
                             is_fvar<scalar_type_t<std::decay_t<T>>>::value> {};

template <typename T>
using enable_if_contains_fvar = std::enable_if_t<is_contains_fvar<T>::value>;

template <typename T>
using enable_if_not_contains_fvar = enable_if_not<is_contains_fvar<T>>;

template <typename... Types>
using enable_if_all_contains_fvar = enable_if_all<is_contains_fvar<Types>...>;

template <typename... Types>
using enable_if_any_contains_fvar = enable_if_any<is_contains_fvar<Types>...>;

template <typename... Types>
using enable_if_all_not_contains_fvar = enable_if_all_not<is_contains_fvar<Types>...>;

template <typename... Types>
using enable_if_any_not_contains_fvar = enable_if_any_not<is_contains_fvar<Types>...>;


// Enables if type is var or arithmetic
template <typename T>
using enable_if_var_or_arithmetic = std::enable_if_t<is_var_or_arithmetic<T>::value>;

template <typename T>
using enable_if_not_var_or_arithmetic = enable_if_not<is_var_or_arithmetic<T>>;

template <typename... Types>
using enable_if_all_var_or_arithmetic = enable_if_all<is_var_or_arithmetic<Types>...>;

template <typename... Types>
using enable_if_any_var_or_arithmetic = enable_if_any<is_var_or_arithmetic<Types>...>;

template <typename... Types>
using enable_if_all_not_var_or_arithmetic = enable_if_all_not<is_var_or_arithmetic<Types>...>;

template <typename... Types>
using enable_if_any_not_var_or_arithmetic = enable_if_any_not<is_var_or_arithmetic<Types>...>;

// Enables if type is var or arithmetic
template <typename T>
struct is_contains_var_or_arithmetic : std::integral_constant<bool, is_var_or_arithmetic<scalar_type_t<T>>::value> {};

template <typename T>
using enable_if_contains_var_or_arithmetic = std::enable_if_t<is_contains_var_or_arithmetic<T>::value>;

template <typename T>
using enable_if_not_contains_var_or_arithmetic = enable_if_not<is_contains_var_or_arithmetic<T>>;

template <typename... Types>
using enable_if_all_contains_var_or_arithmetic = enable_if_all<is_contains_var_or_arithmetic<Types>...>;

template <typename... Types>
using enable_if_any_contains_var_or_arithmetic = enable_if_any<is_contains_var_or_arithmetic<Types>...>;

template <typename... Types>
using enable_if_all_not_contains_var_or_arithmetic = enable_if_all_not<is_contains_var_or_arithmetic<Types>...>;

template <typename... Types>
using enable_if_any_not_contains_var_or_arithmetic = enable_if_any_not<is_contains_var_or_arithmetic<Types>...>;


// Checks whether type is arithmetic, var, or fvar
template <typename T>
struct is_stan_scalar
    : std::integral_constant<bool, std::is_arithmetic<std::decay_t<T>>::value
                                       || is_var<std::decay_t<T>>::value
                                       || is_fvar<std::decay_t<T>>::value> {};

template <typename T>
using enable_if_stan_scalar = std::enable_if_t<is_stan_scalar<T>::value>;

template <typename T>
using enable_if_not_stan_scalar = enable_if_not<is_stan_scalar<T>>;

template <typename... Types>
using enable_if_all_stan_scalar = enable_if_all<is_stan_scalar<Types>...>;

template <typename... Types>
using enable_if_any_stan_scalar = enable_if_any<is_stan_scalar<Types>...>;

template <typename... Types>
using enable_if_all_not_stan_scalar = enable_if_all_not<is_stan_scalar<Types>...>;

template <typename... Types>
using enable_if_any_not_stan_stan_scalar = enable_if_any_not<is_stan_scalar<Types>...>;

// Check whether a type contains (or is) arithmetic, var, or fvar
template <typename T>
struct is_contains_stan_scalar
    : std::integral_constant<
          bool, is_stan_scalar<scalar_type_t<T>>::value> {};

template <typename T>
using enable_if_contains_stan_scalar = std::enable_if_t<is_contains_stan_scalar<T>::value>;

template <typename T>
using enable_if_not_contains_stan_scalar
    = enable_if_not<is_contains_stan_scalar<T>>;

template <typename... Types>
using enable_if_all_contains_stan_scalar = enable_if_all<is_contains_stan_scalar<Types>...>;

template <typename... Types>
using enable_if_any_contains_stan_scalar = enable_if_any<is_contains_stan_scalar<Types>...>;

template <typename... Types>
using enable_if_all_not_contains_stan_scalar = enable_if_all<is_contains_stan_scalar<Types>...>;

template <typename... Types>
using enable_if_any_not_stan_contains_stan_scalar = enable_if_any<is_contains_stan_scalar<Types>...>;




// Checks whether type is a scalar as defined by the standard
template <typename T>
struct is_scalar_decay
    : std::integral_constant<bool, std::is_scalar<std::decay_t<T>>::value> {};

template <typename T>
using enable_if_scalar = std::enable_if_t<is_scalar_decay<T>::value>;

template <typename T>
using enable_if_not_scalar = enable_if_not<is_scalar_decay<T>>;

template <typename... Types>
using enable_if_all_scalar = enable_if_all<is_scalar_decay<Types>...>;

template <typename... Types>
using enable_if_any_scalar = enable_if_any<is_scalar_decay<Types>...>;

template <typename... Types>
using enable_if_all_not_scalar = enable_if_all_not<is_scalar_decay<Types>...>;

template <typename... Types>
using enable_if_any_not_scalar = enable_if_any_not<is_scalar_decay<Types>...>;

}  // namespace stan
#endif
