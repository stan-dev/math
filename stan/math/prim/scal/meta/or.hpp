#ifndef STAN_MATH_PRIM_SCAL_META_OR_HPP
#define STAN_MATH_PRIM_SCAL_META_OR_HPP

#include <type_traits>

namespace stan {
namespace math {

template <typename... Conds>
struct or_ : std::false_type {};

template <typename Cond, typename... Conds>
struct or_<Cond, Conds...>
    : std::conditional<Cond::value, std::true_type, or_<Conds...>>::type {};

template <typename... Conds>
struct or_not_ : std::false_type {};

template <typename Cond, typename... Conds>
struct or_not_<Cond, Conds...>
    : std::conditional<!Cond::value, std::true_type, or_not_<Conds...>>::type {
};

}  // namespace math
}  // namespace stan
#endif
