#ifndef STAN_MATH_PRIM_FUNCTOR_TUPLE_CONCAT_HPP
#define STAN_MATH_PRIM_FUNCTOR_TUPLE_CONCAT_HPP

#include <stan/math/prim/functor/make_holder_tuple.hpp>
#include <stan/math/prim/meta.hpp>
#include <functional>
#include <tuple>
#include <utility>

/**
 * `tuple_concat` only exists because of a bug in clang-7's `std::tuple_cat`
 * If we move up to clang-8+, we can remove these functions and use
 * `std::tuple_cat`
 */

namespace stan {
namespace math {
namespace internal {

template <typename Tuple1, typename Tuple2, std::size_t... I1,
          std::size_t... I2>
inline auto constexpr tuple_concat_impl(Tuple1&& x, Tuple2&& y,
                                        std::index_sequence<I1...> /* i */,
                                        std::index_sequence<I2...> /* j */) {
  return make_holder_tuple(std::get<I1>(std::forward<Tuple1>(x))...,
                           std::get<I2>(std::forward<Tuple2>(y))...);
}

template <typename Tuple1, typename Tuple2, typename Tuple3, std::size_t... I1,
          std::size_t... I2, std::size_t... I3>
inline auto constexpr tuple_concat_impl(Tuple1&& x, Tuple2&& y, Tuple3&& z,
                                        std::index_sequence<I1...> /* i */,
                                        std::index_sequence<I2...> /* j */,
                                        std::index_sequence<I3...> /* k */) {
  return make_holder_tuple(std::get<I1>(std::forward<Tuple1>(x))...,
                           std::get<I2>(std::forward<Tuple2>(y))...,
                           std::get<I3>(std::forward<Tuple3>(z))...);
}
}  // namespace internal

/**
 * Base case to pass a tuple forward.
 * @tparam Tuple Tuple type.
 * @param x Tuple.
 */
inline constexpr auto tuple_concat() noexcept { return std::make_tuple(); }

/**
 * Base case to pass a tuple forward.
 * @tparam Tuple Tuple type.
 * @param x Tuple.
 */
template <typename Tuple>
inline auto tuple_concat(Tuple&& x) noexcept {
  return std::forward<Tuple>(x);
}

/**
 * Concatenates two tuples
 * @tparam Tuple1 First tuple type
 * @tparam Tuple2 Second tuple type
 * @param x First tuple
 * @param y Second tuple
 * @return A tuple containing the elements of x followed by the elements of y
 */
template <typename Tuple1, typename Tuple2>
inline auto tuple_concat(Tuple1&& x, Tuple2&& y) {
  return internal::tuple_concat_impl(
      std::forward<Tuple1>(x), std::forward<Tuple2>(y),
      std::make_index_sequence<std::tuple_size<std::decay_t<Tuple1>>{}>{},
      std::make_index_sequence<std::tuple_size<std::decay_t<Tuple2>>{}>{});
}

/**
 * Concatenates three tuples.
 * @tparam Tuple1 First tuple type
 * @tparam Tuple2 Second tuple type
 * @tparam Tuple3 Third tuple type
 * @param x First tuple
 * @param y Second tuple
 * @param z Third tuple
 * @return A tuple containing the elements of x followed by the elements of y
 * and z
 */
template <typename Tuple1, typename Tuple2, typename Tuple3>
inline auto tuple_concat(Tuple1&& x, Tuple2&& y, Tuple3&& z) {
  return internal::tuple_concat_impl(
      std::forward<Tuple1>(x), std::forward<Tuple2>(y), std::forward<Tuple3>(z),
      std::make_index_sequence<std::tuple_size<std::decay_t<Tuple1>>{}>{},
      std::make_index_sequence<std::tuple_size<std::decay_t<Tuple2>>{}>{},
      std::make_index_sequence<std::tuple_size<std::decay_t<Tuple3>>{}>{});
}

/**
 * Concatenates multiple tuples.
 * @tparam Tuple1 First tuple type
 * @tparam Tuple2 Second tuple type
 * @tparam OtherTuples Remaining tuple types
 * @param x First tuple
 * @param y Second tuple
 * @param args Remaining tuples
 * @return A tuple containing the elements of x followed by the elements of y
 * and the remaining tuples
 */
template <typename Tuple1, typename Tuple2, typename... OtherTuples>
inline auto tuple_concat(Tuple1&& x, Tuple2&& y, OtherTuples&&... args) {
  return tuple_concat(
      tuple_concat(std::forward<Tuple1>(x), std::forward<Tuple2>(y)),
      std::forward<OtherTuples>(args)...);
}

/**
 * Concatenates multiple tuples.
 * @tparam Tuple1 First tuple type
 * @tparam Tuple2 Second tuple type
 * @tparam Tuple3 Third tuple type
 * @tparam OtherTuples Remaining tuple types
 * @param x First tuple
 * @param y Second tuple
 * @param z Third tuple
 * @param args Remaining tuples
 * @return A tuple containing the elements of x followed by the elements of y,
 * z, and the remaining tuples
 */
template <typename Tuple1, typename Tuple2, typename Tuple3,
          typename... OtherTuples>
inline auto tuple_concat(Tuple1&& x, Tuple2&& y, Tuple3&& z,
                         OtherTuples&&... args) {
  return tuple_concat(
      tuple_concat(std::forward<Tuple1>(x), std::forward<Tuple2>(y),
                   std::forward<Tuple3>(z)),
      std::forward<OtherTuples>(args)...);
}
}  // namespace math
}  // namespace stan

#endif
