#ifndef STAN_MATH_PRIM_META_CONTAINS_AUTODIFF_HPP
#define STAN_MATH_PRIM_META_CONTAINS_AUTODIFF_HPP

#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_autodiff.hpp>
#include <tuple>
#include <type_traits>
#include <vector>

namespace stan {

/**
 * Trait indicating whether a type contains at least one autodiff scalar.
 *
 * Scalars are checked directly. Tuple and std::vector specializations recurse
 * into their contained element types.
 *
 * @tparam T type to inspect
 * @return Inherits `bool_constant<true>` when `T` is an autodiff scalar;
 * otherwise `bool_constant<false>`.
 * @throw None.
 */
template <typename T>
struct contains_autodiff : bool_constant<is_autodiff_v<std::decay_t<T>>> {};

/**
 * Tuple specialization. Returns true if any top-level tuple element recursively
 * contains an autodiff scalar.
 *
 * @tparam Args tuple element types
 * @return Inherits `bool_constant<true>` when any tuple element recursively
 * contains an autodiff scalar.
 * @throw None.
 */
template <typename... Args>
struct contains_autodiff<std::tuple<Args...>>
    : bool_constant<(contains_autodiff<std::decay_t<Args>>::value || ...)> {};

/**
 * `std::vector` specialization. Returns whether the element type recursively
 * contains an autodiff scalar.
 *
 * @tparam T vector element type
 * @tparam VecArgs allocator and container argument types
 * @return Inherits `bool_constant<true>` when `T` recursively contains an
 * autodiff scalar.
 * @throw None.
 */
template <typename T, typename... VecArgs>
struct contains_autodiff<std::vector<T, VecArgs...>>
    : bool_constant<contains_autodiff<std::decay_t<T>>::value> {};

/**
 * Convenience variable template for `contains_autodiff`.
 *
 * @tparam T type to inspect
 * @return `true` when `T` recursively contains at least one autodiff scalar.
 * @throw None.
 */
template <typename T>
inline constexpr bool contains_autodiff_v
    = contains_autodiff<std::decay_t<T>>::value;

}  // namespace stan

#endif
