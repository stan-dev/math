#ifndef STAN_MATH_PRIM_FUNCTOR_WALK_TUPLES_HPP
#define STAN_MATH_PRIM_FUNCTOR_WALK_TUPLES_HPP

#include <functional>
#include <tuple>
#include <utility>

namespace stan {
namespace math {
namespace internal {
/**
 * Apply a specified functor to the values from multiple tuples at a
 * provided index
 *
 * @tparam Function Type of functor to apply
 * @tparam Tuples Types of the provided tuples to index
 * @param func Functor to apply
 * @param tuples Tuples to index and apply function to
 */
template <std::size_t Index, typename Function, typename... Tuples>
constexpr void invoke_on_element(Function&& func, Tuples&&... tuples) {
  func(std::get<Index>(std::forward<Tuples>(tuples))...);
}

/**
 * Implementation function for calling invoke_on_element with every index
 * in the provided tuples. The std::initializer trick is used as a pre-C++17
 * alternative for fold expressions.
 *
 * @tparam F Type of provided functor
 * @tparam Tuples Types of tuples to iterate over
 * @param f Functor to apply at each index
 * @param ts Tuples to iterate over
 */
template <class F, std::size_t... I, typename... Tuples>
constexpr void walk_tuples_impl(F&& f, std::index_sequence<I...> i,
                                Tuples&&... ts) {
  (void)std::initializer_list<int>{
      (invoke_on_element<I>(std::forward<F>(f), std::forward<Tuples>(ts)...),
       void(), 0)...};
}
}  // namespace internal

/**
 * Iterate over multiple tuples, using the values at each index as arguments
 * to a provided functor. This function is called for its side-effects, and
 * does not return a value.
 *
 * @tparam F Type of provided functor
 * @tparam Tuples Type of tuples to iterate over
 * @param f Functor to apply
 * @param ts Tuple(s) to use as arguments
 */
template <class F, typename... Tuples>
constexpr void walk_tuples(F&& f, Tuples&&... ts) {
  constexpr auto length
      = std::min({std::tuple_size<std::remove_reference_t<Tuples>>::value...});
  internal::walk_tuples_impl(std::forward<F>(f),
                             std::make_index_sequence<length>{},
                             std::forward<Tuples>(ts)...);
}

}  // namespace math
}  // namespace stan

#endif
