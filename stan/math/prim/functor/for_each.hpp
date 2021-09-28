#ifndef STAN_MATH_PRIM_FUNCTOR_FOR_EACH_HPP
#define STAN_MATH_PRIM_FUNCTOR_FOR_EACH_HPP

#include <stan/math/prim/meta.hpp>
#include <functional>
#include <tuple>
#include <utility>

namespace stan {
namespace math {
namespace internal {

/**
 * Implementation of for_each.
 * @note The static cast to void is used in boost::hana's for_each impl
 *  and is used to suppress unused value warnings from the compiler.
 */
template <typename F, typename T, size_t... Is>
constexpr inline auto for_each(F&& f, T&& t, std::index_sequence<Is...>) {
  using Swallow = int[];
  static_cast<void>(Swallow{
      (static_cast<void>(std::forward<F>(f)(std::get<Is>(std::forward<T>(t)))),
       0)...});
}

/**
 * Implementation of Binary for_each.
 * @note The static cast to void is used in boost::hana's for_each impl
 *  and is used to suppress unused value warnings from the compiler.
 */
template <typename F, typename T1, typename T2, size_t... Is>
constexpr inline auto for_each(F&& f, T1&& t1, T2&& t2,
                               std::index_sequence<Is...>) {
  using Swallow = int[];
  static_cast<void>(Swallow{(
      static_cast<void>(std::forward<F>(f)(std::get<Is>(std::forward<T1>(t1)),
                                           std::get<Is>(std::forward<T2>(t2)))),
      0)...});
}

/**
 * Implementation of ternary for_each.
 * @note The static cast to void is used in boost::hana's for_each impl
 *  and is used to suppress unused value warnings from the compiler.
 */
template <typename F, typename T1, typename T2, typename T3, size_t... Is>
constexpr inline auto for_each(F&& f, T1&& t1, T2&& t2, T3&& t3,
                               std::index_sequence<Is...>) {
  using Swallow = int[];
  static_cast<void>(Swallow{(
      static_cast<void>(std::forward<F>(f)(std::get<Is>(std::forward<T1>(t1)),
                                           std::get<Is>(std::forward<T2>(t2)),
                                           std::get<Is>(std::forward<T3>(t3)))),
      0)...});
}
}  // namespace internal

/**
 * Apply a function to each element of a tuple
 * @tparam F type with a valid `operator()`
 * @tparam T Tuple
 * @param f A functor to apply over each element of the tuple.
 * @param t A tuple
 */
template <typename F, typename T>
constexpr inline auto for_each(F&& f, T&& t) {
  return internal::for_each(
      std::forward<F>(f), std::forward<T>(t),
      std::make_index_sequence<std::tuple_size<std::decay_t<T>>::value>());
}

/**
 * Apply a function to each element of two tuples
 * @tparam F type with a valid `operator()`
 * @tparam T1 Tuple
 * @tparam T2 Another tuple
 * @param f A functor to apply over each element of the tuple.
 * @param t1 A tuple
 * @param t2 Another tuple
 */
template <typename F, typename T1, typename T2>
constexpr inline auto for_each(F&& f, T1&& t1, T2&& t2) {
  constexpr auto t1_size = std::tuple_size<std::decay_t<T1>>::value;
  constexpr auto t2_size = std::tuple_size<std::decay_t<T2>>::value;
  static_assert(t1_size == t2_size,
                "Size Mismatch between t1 and t2 in for_each");
  return internal::for_each(std::forward<F>(f), std::forward<T1>(t1),
                            std::forward<T2>(t2),
                            std::make_index_sequence<t1_size>());
}

/**
 * Apply a function to each element of three tuples
 * @tparam F type with a valid `operator()`
 * @tparam T1 Tuple
 * @tparam T2 Another tuple
 * @tparam T3 Another tuple
 * @param f A functor to apply over each element of the tuple.
 * @param t1 A tuple
 * @param t2 Another tuple
 * @param t3 Another tuple
 */
template <typename F, typename T1, typename T2, typename T3>
constexpr inline auto for_each(F&& f, T1&& t1, T2&& t2, T3&& t3) {
  constexpr auto t1_size = std::tuple_size<std::decay_t<T1>>::value;
  constexpr auto t2_size = std::tuple_size<std::decay_t<T2>>::value;
  constexpr auto t3_size = std::tuple_size<std::decay_t<T3>>::value;
  static_assert(t1_size == t2_size,
                "Size Mismatch between t1 and t2 in for_each");
  static_assert(t1_size == t3_size,
                "Size Mismatch between t1 and t3 in for_each");
  return internal::for_each(std::forward<F>(f), std::forward<T1>(t1),
                            std::forward<T2>(t2), std::forward<T3>(t3),
                            std::make_index_sequence<t1_size>());
}

}  // namespace math
}  // namespace stan

#endif
