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
 */
template <typename T>
struct contains_autodiff : bool_constant<is_autodiff_v<std::decay_t<T>>> {};

template <typename... Args>
struct contains_autodiff<std::tuple<Args...>>
    : bool_constant<(contains_autodiff<std::decay_t<Args>>::value || ...)> {};

template <typename T, typename... VecArgs>
struct contains_autodiff<std::vector<T, VecArgs...>>
    : bool_constant<contains_autodiff<std::decay_t<T>>::value> {};

template <typename T>
inline constexpr bool contains_autodiff_v
    = contains_autodiff<std::decay_t<T>>::value;

}  // namespace stan

#endif
