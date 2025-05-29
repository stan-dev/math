#ifndef STAN_MATH_PRIM_FUNCTOR_MAP_IF_HPP
#define STAN_MATH_PRIM_FUNCTOR_MAP_IF_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/make_holder_tuple.hpp>
#include <stan/math/prim/meta.hpp>
#include <functional>
#include <tuple>
#include <utility>

namespace stan {
namespace math {
namespace internal {

template <template <typename...> class Filter, typename F, typename Arg>
inline constexpr decltype(auto) filter_fun(F&& f, Arg&& arg) {
  if constexpr (Filter<Arg>::value) {
    return f(std::forward<Arg>(arg));
  } else {
    return std::forward<Arg>(arg);
  }
}

}  // namespace internal

/*
 * Subset a tuple by a compile time filter on the types and return a
 * tuple with an unary functor applied to each element where the filter was
 * true:
 * @note For types that fail the filter, a reference to the original type is
 *  returned so we avoid a copy.
 *
 * @tparam Filter a struct that accepts one template parameter and has a static
 *   constexpr bool member named value
 * @tparam F Type of functor
 * @tparam Args A parameter pack of arguments
 * @param f functor callable
 * @param args parameter pack of args
 */
template <template <typename...> class Filter, typename F, typename Tuple,
          require_t<is_tuple<Tuple>>* = nullptr>
inline constexpr auto map_if(F&& f, Tuple&& arg) {
  return stan::math::apply(
      [](auto&& f, auto&&... args) {
        return make_holder_tuple(internal::filter_fun<Filter>(
            std::forward<decltype(f)>(f),
            std::forward<decltype(args)>(args))...);
      },
      std::forward<Tuple>(arg), std::forward<F>(f));
}

/*
 * Subset a parameter pack by a compile time filter on the types and return a
 * tuple with an unary functor applied to each element where the filter was
 * true:
 *
 * @note For types that fail the filter, a reference to the original type is
 *  returned so we avoid a copy.
 * @tparam Filter a struct that accepts one template parameter and has a static
 *   constexpr bool member named value
 * @tparam F Type of functor
 * @tparam Args A parameter pack of arguments
 * @param f functor callable
 * @param args parameter pack of args
 */
template <template <typename...> class Filter, typename F, typename Arg1,
          typename... Args,
          require_t<bool_constant<!is_tuple<Arg1>::value>>* = nullptr>
inline constexpr auto map_if(F&& f, Arg1&& arg1, Args&&... args) {
  return make_holder_tuple(
      internal::filter_fun<Filter>(f, std::forward<Arg1>(arg1)),
      internal::filter_fun<Filter>(f, std::forward<Args>(args))...);
}

template <template <typename...> class Filter, typename F, typename Arg,
          require_t<bool_constant<!is_tuple<Arg>::value>>* = nullptr>
inline constexpr decltype(auto) map_if(F&& f, Arg&& arg) {
  return internal::filter_fun<Filter>(std::forward<F>(f),
                                      std::forward<Arg>(arg));
}

}  // namespace math
}  // namespace stan

#endif
