#ifndef STAN_MATH_PRIM_SCAL_FUN_IS_NAN_HPP
#define STAN_MATH_PRIM_SCAL_FUN_IS_NAN_HPP

#include <utility>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns true if the input is NaN and false otherwise.
 *
 * Delegates to <code>std::isnan</code>.
 *
 * @param x Value to test.
 * @return <code>true</code> if the value is NaN.
 */
template <typename T>
inline bool is_nan(T x) {
  return std::isnan(x);
}

/**
 * Returns true if the input is NaN and false otherwise.
 *
 * Delegates to <code>std::isnan</code>.
 *
 * @param x first argument
 * @param xs parameter pack of remaining arguments to forward to function
 * @return <code>true</code> if any value is NaN
 */
template <typename T, typename... Ts>
inline bool is_nan(T x, Ts... xs) {
  return (is_nan(x) || is_nan(std::forward<Ts>(xs)...));
}
}  // namespace math
}  // namespace stan

#endif
