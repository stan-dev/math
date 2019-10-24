#ifndef STAN_MATH_PRIM_SCAL_FUN_FINITE_DIFF_STEPSIZE_HPP
#define STAN_MATH_PRIM_SCAL_FUN_FINITE_DIFF_STEPSIZE_HPP

#include <cmath>
#include <limits>

namespace stan {
namespace math {

/**
 * Return the stepsize for finite difference evaluations at the
 * specified scalar.
 *
 * <p>The forumula used is `stepsize(u) = cbrt(epsilon) * max(1,
 * abs(u)).`
 *
 * @param u initial value to increment
 * @return stepsize away from u for finite differences
 */
inline double finite_diff_stepsize(double u) {
  static const double cbrt_epsilon
      = std::cbrt(std::numeric_limits<double>::epsilon());
  return cbrt_epsilon * std::fmax(1, fabs(u));
}
}  // namespace math
}  // namespace stan
#endif
