#ifndef STAN_MATH_PRIM_FUN_FINITE_DIFF_STEPSIZE_HPP
#define STAN_MATH_PRIM_FUN_FINITE_DIFF_STEPSIZE_HPP

#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the stepsize for finite difference evaluations at the
 * specified scalar.
 *
 * <p>The formula used is `stepsize(u) = cbrt(epsilon) * max(1,
 * abs(u)).`
 *
 * @param u initial value to increment
 * @return stepsize away from u for finite differences
 */
inline double finite_diff_stepsize(double u) {
  using std::fabs;
  static const double eps = std::numeric_limits<double>::epsilon();
  static const double eps_pow = std::pow(eps, 1.0 / 7.0);   // for 6th-order stencil
  return eps_pow * std::fmax(1.0, fabs(u));
}

}  // namespace math
}  // namespace stan
#endif
