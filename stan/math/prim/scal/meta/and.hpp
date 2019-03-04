#ifndef STAN_MATH_PRIM_SCAL_META_AND_HPP
#define STAN_MATH_PRIM_SCAL_META_AND_HPP

#include <type_traits>

namespace stan {
namespace math {
/**
 * A helper metaprogram to peform conjuction
 * on type conditionals for zero or more provided types.
 * See \link prim/scal/meta/is_var_or_arithmetic.hpp \endlink
 * for an example.
 */
template <typename... T>
struct and_ : std::true_type {};

template <typename T, typename... Ts>
struct and_<T, Ts...>
    : std::conditional<T::value, and_<Ts...>, std::false_type>::type {};

}  // namespace math
}  // namespace stan
#endif
