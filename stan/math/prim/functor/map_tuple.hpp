#ifndef STAN_MATH_PRIM_FUNCTOR_MAP_TUPLE_HPP
#define STAN_MATH_PRIM_FUNCTOR_MAP_TUPLE_HPP

#include <functional>
#include <tuple>
#include <utility>

namespace stan {
namespace math {
namespace internal {


template <class F, typename Tuple, size_t... Is>
constexpr decltype(auto) map_tuple_impl(F&& f, Tuple&& t, std::index_sequence<Is...>) {
  return std::make_tuple(
    f(std::forward<decltype(std::get<Is>(t))>(std::get<Is>(t)))...
  );
}
}  // namespace internal

template <class F, typename TupleT>
constexpr decltype(auto) map_tuple(F&& f, TupleT&& t) {
  return internal::map_tuple_impl(
        std::forward<F>(f),
        std::forward<TupleT>(t),
        std::make_index_sequence<
          std::tuple_size<std::remove_reference_t<TupleT>>{}>{});
}

}  // namespace math
}  // namespace stan

#endif
