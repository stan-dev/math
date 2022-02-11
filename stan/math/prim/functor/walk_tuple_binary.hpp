#ifndef STAN_MATH_PRIM_FUNCTOR_WALK_TUPLE_BINARY_HPP
#define STAN_MATH_PRIM_FUNCTOR_WALK_TUPLE_BINARY_HPP

#include <functional>
#include <tuple>
#include <utility>

namespace stan {
namespace math {
namespace internal {

template <class F, typename Tuple1, typename Tuple2, size_t... Is>
void walk_tuple_binary_impl(F f, Tuple1 t1, Tuple2 t2, std::index_sequence<Is...>) {
  f(std::get<Is>(t1), std::get<Is>(t2))...;
}

}  // namespace internal

template <class F, typename Tuple1, typename Tuple2>
void walk_tuple_binary(F f, const Tuple1& t1, const Tuple2& t2) {
    internal::walk_tuple_binary_impl(
        f, t1, t2, std::make_index_sequence<std::tuple_size<Tuple1>::value>{});
}

}  // namespace math
}  // namespace stan

#endif
