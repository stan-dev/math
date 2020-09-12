#ifndef STAN_MATH_PRIM_FUN_IS_ANY_NAN_HPP
#define STAN_MATH_PRIM_FUN_IS_ANY_NAN_HPP

#include <stan/math/prim/fun/is_nan.hpp>
#include <utility>

namespace stan {
namespace math {

/**
 * Returns true if the input is NaN and false otherwise.
 *
 * Delegates to <code>stan::math::is_nan</code> so that
 * appropriate specializations can be loaded for autodiff
 * types.
 *
 * @param x Value to test.
 * @return <code>true</code> if the value is NaN.
 */
template <typename T>
inline bool is_any_nan(const T& x) {
  return is_nan(x);
}

/**
 * Returns <code>true</code> if any input is NaN and false otherwise.
 *
 * Delegates to <code>stan::math::is_nan</code>.
 *
 * @param x first argument
 * @param xs parameter pack of remaining arguments to forward to function
 * @return <code>true</code> if any value is NaN
 */
template <typename T, typename... Ts>
inline bool is_any_nan(const T& x, const Ts&... xs) {
  return is_any_nan(x) || is_any_nan(xs...);
}

}  // namespace math
}  // namespace stan

#endif
