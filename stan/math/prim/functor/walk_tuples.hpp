#ifndef STAN_MATH_PRIM_FUNCTOR_WALK_TUPLES_HPP
#define STAN_MATH_PRIM_FUNCTOR_WALK_TUPLES_HPP

#include <functional>
#include <tuple>
#include <utility>

namespace stan {
namespace math {
namespace internal {

template <std::size_t Index, typename Function, typename... Tuples>
constexpr void call_elem(Function&& func, Tuples&&... tuples) {
    func(std::get<Index>(std::forward<Tuples>(tuples))...);
}

template <class F, std::size_t... I, typename... Tuples>
void walk_tuples_impl(F f, std::index_sequence<I...> i,
                      const Tuples&... ts) {
  std::initializer_list<int>{ (call_elem<I>(f, ts...), 0)... };
}

}  // namespace internal

template <class F, typename... Tuples>
void walk_tuples(F f, const Tuples&... ts) {
  constexpr auto length = std::min({std::tuple_size<std::remove_reference_t<Tuples>>::value...});
  internal::walk_tuples_impl(f, std::make_index_sequence<length>{}, ts...);
}


}  // namespace math
}  // namespace stan

#endif
