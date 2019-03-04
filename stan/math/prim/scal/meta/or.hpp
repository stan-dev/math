#ifndef STAN_MATH_PRIM_SCAL_META_OR_HPP
#define STAN_MATH_PRIM_SCAL_META_OR_HPP

#include <type_traits>

namespace stan {

/**
 * A helper metaprogram to peform disjunction
 * on type conditionals for zero or more provided types.
 * See \link prim/scal/meta/contains_vector.hpp \endlink
 * for an example.
 */
template <typename... Conds>
struct or_ : std::false_type {};

template <typename Cond, typename... Conds>
struct or_<Cond, Conds...>
    : std::conditional<Cond::value, std::true_type, or_<Conds...>>::type {};
/**
 * A helper metaprogram to peform disjunction
 * on negated type conditionals for zero or
 * more provided types.
 * See \link prim/scal/meta/contains_vector.hpp \endlink
 * for an example.
 */
template <typename... Conds>
struct or_not_ : std::false_type {};

template <typename Cond, typename... Conds>
struct or_not_<Cond, Conds...>
    : std::conditional<!Cond::value, std::true_type, or_not_<Conds...>>::type {
};

}  // namespace stan
#endif
