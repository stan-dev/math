#ifndef STAN_MATH_PRIM_SCAL_FUN_LGAMMA_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LGAMMA_HPP

#include <cmath>

namespace stan {
namespace math {

/**
 * Return the natural logarithm of the gamma function applied to the
 * specified argument.  If the input is a non-positive integer, the
 * result is positive infinity.  If the input is not-a-number, the
 * result and gradient are not-a-number.
 *
 * @param x argument
 * @return natural logarithm of the gamma function applied to
 * argument
 */
inline double lgamma(double x) { return std::lgamma(x); }

/**
 * Return the natural logarithm of the gamma function applied
 * to the specified argument.
 *
 * @param x argument
 * @return natural logarithm of the gamma function applied to
 * argument
 */
inline double lgamma(int x) { return std::lgamma(x); }

}  // namespace math
}  // namespace stan
#endif
