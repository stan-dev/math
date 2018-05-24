#ifndef STAN_MATH_PRIM_SCAL_FUN_STEP_HPP
#define STAN_MATH_PRIM_SCAL_FUN_STEP_HPP

namespace stan {
namespace math {

/**
 * The step, or Heaviside, function.
 *
 * For double NaN input, step(NaN) returns 0.
 *
 * Note: behavior changed from
 * <code>step(y) = (y < 0.0) ? 0 : 1</code>
 * to
 *
   \f[
   \mbox{step}(y) =
   \begin{cases}
     0 & \mbox{if } y \leq 0 \\
     1 & \mbox{if } y > 0  \\[6pt]
     0 & \mbox{if } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam T value type
 * @param[in] y value
 * @return 1 if value is greater than 0 and 0 otherwise
 */
template <typename T>
inline double step(const T& y) {
  return y > 0.0;
}

}  // namespace math
}  // namespace stan

#endif
