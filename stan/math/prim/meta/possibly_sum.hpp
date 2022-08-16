#ifndef STAN_MATH_PRIM_META_POSSIBLY_SUM_HPP
#define STAN_MATH_PRIM_META_POSSIBLY_SUM_HPP

#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>

namespace stan {
namespace math {

/**
 * Conditionally sum the input at compile time.
 * @tparam CondSum A struct with a static boolean member `value` which if true
 *  will allow the input value to be summed
 * @tparam T A scalar, Eigen type, or standard vector with inner scalar type.
 * @param x The value to be summed.
 */
template <typename CondSum, typename T, require_t<CondSum>* = nullptr>
inline auto possibly_sum(T&& x) {
  return sum(std::forward<T>(x));
}

/**
 * Conditionally sum the input at compile time. This overload does not sum.
 * @tparam CondSum A struct with a static boolean member `value` which if false
 *  will pass the input to the output.
 * @tparam T A scalar, Eigen type, or standard vector with inner scalar type.
 * @param x The value to be passed trhough.
 */
template <typename CondSum, typename T1, require_not_t<CondSum>* = nullptr>
inline auto possibly_sum(T1&& x) {
  return std::forward<T1>(x);
}

}  // namespace math
}  // namespace stan

#endif
