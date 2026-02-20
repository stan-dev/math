#ifndef STAN_MATH_PRIM_FUNCTOR_FILTER_TYPES_HPP
#define STAN_MATH_PRIM_FUNCTOR_FILTER_TYPES_HPP

#include <stan/math/prim/functor/filter_map.hpp>
#include <stan/math/prim/meta.hpp>
#include <utility>

namespace stan {
namespace math {

/**
 * Filter a tuple and return a tuple containing references to entries whose
 * types satisfy `Predicate`.
 *
 * This helper preserves nested tuple structure and delegates recursion to
 * `filter_map`.
 *
 * @tparam Predicate unary type-trait predicate
 * @tparam T tuple-like input type
 * @param[in] t input tuple to filter
 * @return tuple containing only filtered entries; returned elements preserve
 * reference/value category from `t`
 * @throw None.
 */
template <template <typename...> class Predicate, typename T,
          require_tuple_t<T>* = nullptr>
inline constexpr decltype(auto) filter_types(T&& t) {
  return stan::math::filter_map<Predicate>(
      [](auto&& arg) -> decltype(auto) {
        return std::forward<decltype(arg)>(arg);
      },
      std::forward<T>(t));
}

}  // namespace math
}  // namespace stan

#endif
