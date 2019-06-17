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
  using std::cbrt;
  using std::fmax;
  using std::numeric_limits;
  static const double cbrt_epsilon = cbrt(numeric_limits<double>::epsilon());
  return cbrt_epsilon * fmax(1, fabs(u));
}
}  // namespace math
}  // namespace stan
#endif
