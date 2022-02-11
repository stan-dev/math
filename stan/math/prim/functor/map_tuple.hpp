#ifndef STAN_MATH_PRIM_FUNCTOR_MAP_TUPLE_HPP
#define STAN_MATH_PRIM_FUNCTOR_MAP_TUPLE_HPP

#include <functional>
#include <tuple>
#include <utility>

namespace stan {
namespace math {
namespace internal {


  template <class F, typename Tuple, size_t... Is>
  auto map_tuple_impl(F f, Tuple t, std::index_sequence<Is...>) {
      return std::make_tuple(
          f(std::get<Is>(t))...
      );
  }
}  // namespace internal

  template <class F, typename... Args>
  auto map_tuple(F f, const std::tuple<Args...>& t) {
      return internal::map_tuple_impl(
          f, t, std::make_index_sequence<sizeof...(Args)>{});
  }

}  // namespace math
}  // namespace stan

#endif
