#ifndef STAN_MATH_PRIM_META_INDEX_APPLY_HPP
#define STAN_MATH_PRIM_META_INDEX_APPLY_HPP
#include <utility>
#include <cstddef>

namespace stan {
namespace math {

namespace internal {
template <class F, std::size_t... Is>
constexpr auto index_apply_impl(F f, std::index_sequence<Is...>) {
  return f(std::integral_constant<std::size_t, Is>{}...);
}
}  // namespace internal

/**
 * Calls given callable with an index sequence. Combination of this and a lambda
 * can be used to avoid separating implementation that needs an index sequence
 * from functions that operate on tuples of parameter packs.
 * @tparam N size of index sequence to use
 * @tparam F callable type
 * @param f callable
 * @return what the given callable returns
 */
template <std::size_t N, class F>
constexpr auto index_apply(F f) {
  return internal::index_apply_impl(f, std::make_index_sequence<N>{});
}

}  // namespace math
}  // namespace stan

#endif  // STAN_MATH_PRIM_META_APPLY_HPP
