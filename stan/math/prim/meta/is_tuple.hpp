#ifndef STAN_MATH_PRIM_META_IS_TUPLE_HPP
#define STAN_MATH_PRIM_META_IS_TUPLE_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <tuple>
#include <type_traits>

namespace stan {
namespace math {
namespace internal {

template <typename... Types>
struct is_tuple_impl : std::false_type {};

template <typename... Types>
struct is_tuple_impl<std::tuple<Types...>> : std::true_type {};
}  // namespace internal

template <typename T>
struct is_tuple : internal::is_tuple_impl<std::decay_t<T>> {};

STAN_ADD_REQUIRE_UNARY(tuple, is_tuple, require_std);
STAN_ADD_REQUIRE_UNARY_INNER(tuple, is_tuple, require_std);

}  // namespace math
}  // namespace stan

#endif
