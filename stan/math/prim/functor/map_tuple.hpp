#ifndef STAN_MATH_PRIM_FUNCTOR_MAP_TUPLE_HPP
#define STAN_MATH_PRIM_FUNCTOR_MAP_TUPLE_HPP

#include <functional>
#include <tuple>
#include <utility>

namespace stan {
namespace math {
namespace internal {

/**
 * Invoke the functor f over each argument given in t and indexed in the index
 * sequence I
 *
 * @tparam F Type of functor
 * @tparam Tuple Type of tuple containing arguments
 * @tparam I Parameter pack of an index sequence going from 0 to
 * std::tuple_size<T>::value - 1 inclusive
 * @param f functor callable
 * @param t tuple of arguments
 * @param i placeholder variable for index sequence
 */
template <class F, typename Tuple, size_t... Is>
constexpr decltype(auto) map_tuple_impl(
  F&& f, Tuple&& t, std::index_sequence<Is...>) {
  return std::make_tuple(
    f(std::forward<decltype(std::get<Is>(t))>(std::get<Is>(t)))...);
}
}  // namespace internal

/**
 * Call the functor f over each element of the tuple of arguments t,
 * returning a tuple of the same size:
 *
 * std::make_tuple(f(std::get<0>(t)), f(std::get<1>(t))...)
 *
 * @tparam F Type of functor
 * @tparam Tuple Type of tuple containing arguments
 * @param f functor callable
 * @param t tuple of arguments
 */
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
