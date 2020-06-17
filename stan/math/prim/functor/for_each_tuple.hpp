#ifndef STAN_MATH_PRIM_FUNCTOR_FOR_EACH_TUPLE_HPP
#define STAN_MATH_PRIM_FUNCTOR_FOR_EACH_TUPLE_HPP

#include <stan/math/prim/meta.hpp>
#include <functional>
#include <tuple>
#include <utility>

namespace stan {
namespace math {
namespace internal {
template<typename F, typename T, size_t... Is>
static inline auto for_each(F&& f, T&& t, std::index_sequence<Is...>) {
  auto l = { (std::forward<F>(f)(std::get<Is>(t)), 0)...};
}
template<typename F, typename T1, typename T2, size_t... Is>
static inline auto for_each(F&& f, T1&& t1, T2&& t2, std::index_sequence<Is...>) {
    auto l = { (std::forward<F>(f)(std::get<Is>(t1), std::get<Is>(t2)), 0)...};
}
}
template<typename F, typename T>
static inline auto for_each_tuple(F&& f, T&& t) {
  return internal::for_each(std::forward<F>(f), std::forward<T>(t),
     std::make_index_sequence<std::tuple_size<std::decay_t<T>>::value>());
}
template<typename F, typename T1, typename T2>
static inline auto for_each_tuple(F&& f, T1&& t1, T2&& t2) {
  return internal::for_each(std::forward<F>(f), std::forward<T1>(t1), std::forward<T2>(t2),
     std::make_index_sequence<std::tuple_size<std::decay_t<T1>>::value>());
}

}
}
#endif
