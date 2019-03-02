#ifndef STAN_MATH_PRIM_SCAL_META_AND_HPP
#define STAN_MATH_PRIM_SCAL_META_AND_HPP

#include <type_traits>

namespace stan {
namespace math {

template <typename... Conds>
struct and_ : std::true_type {};

template <typename Cond, typename... Conds>
struct and_<Cond, Conds...>
    : std::conditional<Cond::value, and_<Conds...>, std::false_type>::type {};

}  // namespace math
}  // namespace stan
#endif
