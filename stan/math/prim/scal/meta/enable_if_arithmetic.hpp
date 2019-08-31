#ifndef STAN_MATH_PRIM_SCAL_META_REQUIRE_GENERICS_HPP
#define STAN_MATH_PRIM_SCAL_META_REQUIRE_GENERICS_HPP

#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>

#include <type_traits>

namespace stan {

/**
 * Template metaprogram for deducing if a method should be remove a
 * a function from overload resolution. This method performs a negation of
 * the input boolean.
 */
template <typename Check>
using require_base = std::enable_if_t<Check::value>;

/**
 * Template metaprogram for deducing if a method should be remove a
 * a function from overload resolution. This method performs a negation of
 * the input boolean.
 */
template <typename Check>
using require_not = std::enable_if_t<!Check::value>;

/**
 * Template metaprogram for deducing if a method should be remove a
 * a function from overload resolution.
 *
 * Returns a type void if all conditions are true and otherwise fails.
 */
template <class... Checks>
using require_all = std::enable_if_t<math::conjunction<Checks...>::value>;

/**
 * Template metaprogram for deducing if a method should be remove a
 * a function from overload resolution.
 *
 * Returns a type void if any of the conditions are true and otherwise fails.
 */
template <class... Checks>
using require_any = std::enable_if_t<math::disjunction<Checks...>::value>;

/**
 * Template metaprogram for deducing if a method should be remove a
 * a function from overload resolution.
 *
 * Returns a type void if all of the conditions are false.
 */
template <class... Checks>
using require_not_all = std::enable_if_t<!math::conjunction<Checks...>::value>;

/**
 * Template metaprogram for deducing if a method should be remove a
 * a function from overload resolution.
 *
 * Returns a type void if any of the conditions are false.
 */
template <class... Checks>
using require_not_any = std::enable_if_t<!math::disjunction<Checks...>::value>;

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
using require_not_all_arithmetic
    = require_not_all<std::is_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_not_any_arithmetic
    = require_not_any<std::is_arithmetic<std::decay_t<Types>>...>;

}  // namespace stan
#endif
