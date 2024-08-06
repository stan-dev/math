#ifndef STAN_MATH_PRIM_FUNCTOR_FILTER_MAP_HPP
#define STAN_MATH_PRIM_FUNCTOR_FILTER_MAP_HPP

#include <functional>
#include <tuple>
#include <utility>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/meta.hpp>
namespace stan {
namespace math {
namespace internal {
template <template <typename> class Filter, typename F, typename Arg>
inline constexpr auto filter_fun(F&& f, Arg&& arg) {
  if constexpr (Filter<Arg>::value) {
    return f(std::forward<Arg>(arg));
  } else {
    return std::forward<Arg>(arg);
  }
}
}

/*
 * Subset a tuple by a compile time filter on the types and return a
 * tuple with an unary functor applied to each element where the filter was
 * true:
 *
 * @tparam Filter a struct that accepts one template parameter and has a static
 *   constexpr bool member named value
 * @tparam F Type of functor
 * @tparam Args A parameter pack of arguments
 * @param f functor callable
 * @param args parameter pack of args
 */
template <template <typename> class Filter, typename F, typename Tuple,
          require_t<is_tuple<Tuple>>* = nullptr>
inline constexpr auto filter_map(F&& f, Tuple&& arg) {
  return stan::math::apply([](auto&& f, auto&&... args) {
    return std::make_tuple(
      internal::filter_fun<Filter>(std::forward<decltype(f)>(f), args)...);
  }, std::forward<Tuple>(arg), f);
}

/*
 * Subset a parameter pack by a compile time filter on the types and return a
 * tuple with an unary functor applied to each element where the filter was
 * true:
 *
 * @tparam Filter a struct that accepts one template parameter and has a static
 *   constexpr bool member named value
 * @tparam F Type of functor
 * @tparam Args A parameter pack of arguments
 * @param f functor callable
 * @param args parameter pack of args
 */
template <template <typename> class Filter, typename F, typename Arg1,
          typename... Args,
          require_t<bool_constant<!is_tuple<Arg1>::value>>* = nullptr>
inline constexpr auto filter_map(F&& f, Arg1&& arg1, Args&&... args) {
  return filter_map<Filter>(
    std::forward<F>(f), std::forward_as_tuple(
      std::forward<Arg1>(arg1), std::forward<Args>(args)...));
}

template <template <typename> class Filter, typename F, typename Arg,
          require_t<bool_constant<!is_tuple<Arg>::value>>* = nullptr>
inline constexpr auto filter_map(F&& f, Arg&& arg1) {
  return internal::filter_fun<Filter>(std::forward<F>(f), std::forward<Arg>(arg));
}

}  // namespace math
}  // namespace stan

#endif
