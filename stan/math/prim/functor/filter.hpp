#ifndef STAN_MATH_PRIM_FUNCTOR_FILTER_HPP
#define STAN_MATH_PRIM_FUNCTOR_FILTER_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/partially_forward_as_tuple.hpp>
#include <stan/math/prim/functor/tuple_concat.hpp>
#include <stan/math/prim/meta.hpp>
#include <functional>
#include <tuple>
#include <utility>

namespace stan {
namespace math {

template <template <typename...> class Filter, std::size_t Index = 0,
          typename F, typename Tuple>
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

}  // namespace math
}  // namespace stan

#endif
