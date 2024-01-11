#ifndef STAN_MATH_PRIM_FUNCTOR_APPLY_HPP
#define STAN_MATH_PRIM_FUNCTOR_APPLY_HPP

#include <functional>
#include <tuple>
#include <utility>

namespace stan {
namespace math {
namespace internal {
/*
 * Invoke the functor f with arguments given in t and indexed in the index
 * sequence I with other arguments possibly before or after
 *
 * @tparam F Type of functor
 * @tparam Tuple Type of tuple containing arguments
 * @tparam PreArgs Parameter pack of arguments before the tuple
 * @tparam I Parameter pack of an index sequence going from 0 to
 * std::tuple_size<T>::value - 1 inclusive
 * @param f functor callable
 * @param t tuple of arguments
 * @param i placeholder variable for index sequence
 * @param pre_args parameter pack of arguments to place before elements in
 * tuple.
 */
template <class F, class Tuple, typename... PreArgs, std::size_t... I>
constexpr decltype(auto) apply_impl(F&& f, Tuple&& t,
                                    std::index_sequence<I...> i,
                                    PreArgs&&... pre_args) {
  return std::forward<F>(f)(
      std::forward<PreArgs>(pre_args)...,
      std::forward<decltype(std::get<I>(t))>(std::get<I>(t))...);
}
}  // namespace internal

/*
 * Call the functor f with the tuple of arguments t, like:
 *
 * f(std::get<0>(t), std::get<1>(t), ...)
 *
 * TODO: replace this with implementation in C++ std when C++17 is available
 *
 * @tparam F Type of functor
 * @tparam Tuple Type of tuple containing arguments
 * @tparam PreArgs Parameter pack of arguments before the tuple
 * @param f functor callable
 * @param t tuple of arguments
 * @param pre_args parameter pack of arguments to place before elements in
 * tuple.
 */
template <class F, class Tuple, typename... PreArgs>
constexpr decltype(auto) apply(F&& f, Tuple&& t, PreArgs&&... pre_args) {
  return internal::apply_impl(
      std::forward<F>(f), std::forward<Tuple>(t),
      std::make_index_sequence<
          std::tuple_size<std::remove_reference_t<Tuple>>{}>{},
      std::forward<PreArgs>(pre_args)...);
}

}  // namespace math
}  // namespace stan

#endif
