#ifndef STAN_MATH_PRIM_FUNCTOR_FILTER_MAP_HPP
#define STAN_MATH_PRIM_FUNCTOR_FILTER_MAP_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/partially_forward_as_tuple.hpp>
#include <stan/math/prim/functor/tuple_concat.hpp>
#include <stan/math/prim/meta.hpp>
#include <functional>
#include <tuple>
#include <utility>

namespace stan {
namespace math {

template <template <typename...> class Filter,
          typename F>
inline constexpr auto filter_map(F&& f) noexcept {
  return std::tuple<>{};
}
/**
 * Filter a tuple and apply a functor to each element that passes the filter.
 * @tparam Filter a struct that accepts one template parameter and has a static
 *  constexpr bool member named value that is true if the type should be
 *  included in the output tuple.
 * @tparam Index Index of the current element in the tuple.
 * @tparam F Type of functor
 * @tparam Tuple A tuple
 * @param f functor callable
 * @param tup tuple of arguments
 * @return a tuple with the functor applied to each element which passed the
 * filter.
 */
template <template <typename...> class Filter, std::size_t Index = 0,
          typename F, typename Tuple, require_tuple_t<Tuple>* = nullptr>
inline constexpr auto filter_map(F&& f, Tuple&& tup) {
  if constexpr (Index == (std::tuple_size<std::decay_t<Tuple>>::value)) {
    return std::make_tuple();
  } else {
    constexpr bool apply_filter_b
        = Filter<std::decay_t<decltype(std::get<Index>(tup))>>::value;
    if constexpr (apply_filter_b) {
      if constexpr (stan::math::is_tuple_v<
                        std::tuple_element_t<Index, std::decay_t<Tuple>>>) {
        // This will look like tuple(tuple(tuple(1, 2)), tuple(3, 4)) ->
        // tuple(tuple(1, 2), 3, 4)
        return tuple_concat(std::make_tuple(filter_map<Filter>(
                                f, std::get<Index>(std::forward<Tuple>(tup)))),
                            filter_map<Filter, Index + 1>(
                                std::forward<F>(f), std::forward<Tuple>(tup)));
      } else {
        return tuple_concat(partially_forward_as_tuple(
                                f(std::get<Index>(std::forward<Tuple>(tup)))),
                            filter_map<Filter, Index + 1>(
                                std::forward<F>(f), std::forward<Tuple>(tup)));
      }
    } else {
      return filter_map<Filter, Index + 1>(std::forward<F>(f),
                                           std::forward<Tuple>(tup));
    }
  }
}

template <template <typename...> class Filter,
          typename F, typename T1, typename... Types, require_not_tuple_t<T1>* = nullptr>
inline constexpr auto filter_map(F&& f, T1&& x, Types&&... xs) {
  constexpr bool apply_filter_b = Filter<std::decay_t<T1>>::value;
  if constexpr (apply_filter_b) {
    if (sizeof...(Types) == 0) {
      return partially_forward_as_tuple(std::forward<F>(f)(std::forward<T1>(x)));
    }
    if constexpr (stan::math::is_tuple_v<T1>) {
      // This will look like tuple(tuple(tuple(1, 2)), tuple(3, 4)) ->
      // tuple(tuple(1, 2), 3, 4)
      return tuple_concat(std::make_tuple(filter_map<Filter>(f, std::forward<T1>(x))),
                          filter_map<Filter>(std::forward<F>(f), std::forward<Types>(xs)...));
    } else {
      return tuple_concat(partially_forward_as_tuple(
                              f(std::forward<T1>(x))),
                          filter_map<Filter>(
                              std::forward<F>(f), std::forward<Types>(xs)...));
    }
  } else {
    return filter_map<Filter>(std::forward<F>(f), std::forward<Types>(xs)...);
  }
}


}  // namespace math
}  // namespace stan

#endif