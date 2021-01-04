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
  static_cast<void>(Swallow{(static_cast<void>(std::forward<F>(f)(
                                 std::get<Is>(std::forward<T>(t)), Is)),
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
  static_cast<void>(Swallow{(static_cast<void>(std::forward<F>(f)(
                                 std::get<Is>(std::forward<T1>(t1)),
                                 std::get<Is>(std::forward<T2>(t2)), Is)),
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
  check_size_match("for_each", "t1", std::tuple_size<std::decay_t<T1>>::value,
                   "t2", std::tuple_size<std::decay_t<T2>>::value);
  return internal::for_each(
      std::forward<F>(f), std::forward<T1>(t1), std::forward<T2>(t2),
      std::make_index_sequence<std::tuple_size<std::decay_t<T1>>::value>());
}

}  // namespace math
}  // namespace stan

#endif
