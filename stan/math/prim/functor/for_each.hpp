#ifndef STAN_MATH_PRIM_FUNCTOR_for_each_HPP
#define STAN_MATH_PRIM_FUNCTOR_for_each_HPP

#include <stan/math/prim/meta.hpp>
#include <functional>
#include <tuple>
#include <utility>

namespace stan {
namespace math {
namespace internal {
template<typename F, typename T, size_t... Is>
constexpr inline auto for_each(F&& f, T&& t, std::index_sequence<Is...>) {
  std::initializer_list<int>{(
    static_cast<void>(std::forward<F>(f)(std::get<Is>(std::forward<T>(t)))), 0)...};
}
template<typename F, typename T1, typename T2, size_t... Is>
constexpr inline auto for_each(F&& f, T1&& t1, T2&& t2, std::index_sequence<Is...>) {
    std::initializer_list<int>{(static_cast<void>(std::forward<F>(f)(std::get<Is>(std::forward<T1>(t1)), std::get<Is>(std::forward<T2>(t2)))), 0)...};
}
}
/**
 * Apply a function to each element of a tuple
 * @tparam F type with a valid `operator()`
 * @tparam T tuple
 * @param f A functor to apply over each element of the tuple.
 * @param t A tuple.
 */
template<typename F, typename T>
constexpr inline auto for_each(F&& f, T&& t) {
  return internal::for_each(std::forward<F>(f), std::forward<T>(t),
     std::make_index_sequence<std::tuple_size<std::decay_t<T>>::value>());
}
/**
 * Apply a function to each element of a tuple and other container
 * @tparam F type with a valid `operator()`
 * @tparam T1 tuple
 * @tparam Any container with an `operator[]`
 * @param f A functor to apply over each element of the tuple.
 * @param t A tuple.
 */
template<typename F, typename T1, typename T2>
constexpr inline auto for_each(F&& f, T1&& t1, T2&& t2) {
  return internal::for_each(std::forward<F>(f), std::forward<T1>(t1), std::forward<T2>(t2),
     std::make_index_sequence<std::tuple_size<std::decay_t<T1>>::value>());
}

}
}
#endif
