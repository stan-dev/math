#ifndef STAN_MATH_PRIM_META_IS_DETECTED_HPP
#define STAN_MATH_PRIM_META_IS_DETECTED_HPP

#include <stan/math/prim/meta/void_t.hpp>
#include <type_traits>

namespace stan {

template <typename, template <typename...> class, typename = void>
struct is_detected : std::false_type {};

/**
 * Checks whether a valid type is detected. Most commonly used to detect
 *  attributes of objects.
 * @tparam T The type to be checked.
 * @tparam Op a template template type which attempts to define a member
 *  of an object.
 */
template <typename T, template <typename...> class Op>
struct is_detected<T, Op, void_t<Op<T>>> : std::true_type {};

}  // namespace stan
#endif
