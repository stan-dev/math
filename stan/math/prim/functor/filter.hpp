#ifndef STAN_MATH_PRIM_FUNCTOR_FILTER_HPP
#define STAN_MATH_PRIM_FUNCTOR_FILTER_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/meta.hpp>
#include <functional>
#include <tuple>
#include <utility>

namespace stan {
namespace math {

template <template <typename> class Filter, std::size_t Index = 0,typename F, typename Tuple>
inline constexpr auto filter(F&& f, Tuple&& tup) {
  if constexpr (Index == std::tuple_size<std::decay_t<Tuple>>::value) {
    return std::make_tuple();
  } else if constexpr (Filter<std::tuple_element_t<Index, std::decay_t<Tuple>>>::value) {
    return std::tuple_cat(
      std::tuple(f(std::get<Index>(std::forward<Tuple>(tup)))),
      filter<Filter, Index + 1>(std::forward<F>(f), std::forward<Tuple>(tup)));
  } else {
    return filter<Filter, Index + 1>(std::forward<F>(f), std::forward<Tuple>(tup));
  }
}


}  // namespace math
}  // namespace stan

#endif
