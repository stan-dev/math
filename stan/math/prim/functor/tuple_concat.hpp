#ifndef STAN_MATH_PRIM_FUNCTOR_TUPLE_CONCAT_HPP
#define STAN_MATH_PRIM_FUNCTOR_TUPLE_CONCAT_HPP

#include <stan/math/prim/functor/partially_forward_as_tuple.hpp>
#include <stan/math/prim/meta.hpp>
#include <functional>
#include <tuple>
#include <utility>


namespace stan {
namespace math {
namespace internal {

template <typename Tuple1, typename Tuple2, std::size_t... I1,
          std::size_t... I2>
inline auto constexpr tuple_concat_impl(Tuple1&& x, Tuple2&& y,
                                   std::index_sequence<I1...> i,
                                   std::index_sequence<I2...> j) {
  return std::make_tuple(
      (std::forward<decltype(std::get<I1>(std::forward<Tuple1>(x)))>(
          std::get<I1>(std::forward<Tuple1>(x))))...,
      (std::forward<decltype(std::get<I2>(std::forward<Tuple2>(y)))>(
          std::get<I2>(std::forward<Tuple2>(y))))...);
}

template <typename Tuple1, typename Tuple2, typename Tuple3, std::size_t... I1,
          std::size_t... I2, std::size_t... I3>
inline auto constexpr tuple_concat_impl(Tuple1&& x, Tuple2&& y, Tuple3&& z,
                                   std::index_sequence<I1...> i,
                                   std::index_sequence<I2...> j,
                                   std::index_sequence<I3...> k) {
  return std::make_tuple(
      (std::forward<decltype(std::get<I1>(std::forward<Tuple1>(x)))>(
          std::get<I1>(std::forward<Tuple1>(x))))...,
      (std::forward<decltype(std::get<I2>(std::forward<Tuple2>(y)))>(
          std::get<I2>(std::forward<Tuple2>(y))))...,
      (std::forward<decltype(std::get<I3>(std::forward<Tuple3>(z)))>(
          std::get<I3>(std::forward<Tuple3>(z))))...);
}
}

template <typename Tuple>
inline auto tuple_concat(Tuple&& x) {
  return std::forward<Tuple>(x);
}

template <typename Tuple1, typename Tuple2>
inline auto tuple_concat(Tuple1&& x, Tuple2&& y) {
  return internal::tuple_concat_impl(
      std::forward<Tuple1>(x), std::forward<Tuple2>(y),
      std::make_index_sequence<
          std::tuple_size<std::remove_reference_t<Tuple1>>{}>{},
      std::make_index_sequence<
          std::tuple_size<std::remove_reference_t<Tuple2>>{}>{});
}

template <typename Tuple1, typename Tuple2, typename Tuple3>
inline auto tuple_concat(Tuple1&& x, Tuple2&& y, Tuple3&& z) {
  return internal::tuple_concat_impl(
      std::forward<Tuple1>(x), std::forward<Tuple2>(y), std::forward<Tuple3>(z),
      std::make_index_sequence<
          std::tuple_size<std::remove_reference_t<Tuple1>>{}>{},
      std::make_index_sequence<
          std::tuple_size<std::remove_reference_t<Tuple2>>{}>{},
      std::make_index_sequence<
          std::tuple_size<std::remove_reference_t<Tuple2>>{}>{});
}

template <typename Tuple1, typename Tuple2, typename... OtherTuples>
inline auto tuple_concat(Tuple1&& x, Tuple2&& y, OtherTuples&&... args) {
  return tuple_concat(tuple_concat(std::forward<Tuple1>(x), std::forward<Tuple2>(y)),
                 std::forward<OtherTuples>(args)...);
}

template <typename Tuple1, typename Tuple2, typename Tuple3,
          typename... OtherTuples>
inline auto tuple_concat(Tuple1&& x, Tuple2&& y, Tuple3&& z, OtherTuples&&... args) {
  return tuple_concat(tuple_concat(std::forward<Tuple1>(x), std::forward<Tuple2>(y),
                         std::forward<Tuple3>(z)),
                 std::forward<OtherTuples>(args)...);
}
}
}

#endif
