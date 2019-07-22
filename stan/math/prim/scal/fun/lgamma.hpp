#ifndef STAN_MATH_PRIM_SCAL_FUN_LGAMMA_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LGAMMA_HPP

#if !__MINGW32__
// _REENTRANT must be defined during compilation to ensure that cmath
// exports reentrant safe lgamma_r version.
#include <cmath>
#else
// MinGW does not provide the reentrant lgamma_r such that we fall
// back to boost.
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/boost_policy.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <limits>
#endif

namespace stan {
namespace math {

/**
 * Return the natural logarithm of the gamma function applied to
 * the specified argument.
 *
   \f[
   \mbox{lgamma}(x) =
   \begin{cases}
     \textrm{error} & \mbox{if } x\in \{\dots, -3, -2, -1, 0\}\\
     \ln\Gamma(x) & \mbox{if } x\not\in \{\dots, -3, -2, -1, 0\}\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{lgamma}(x)}{\partial x} =
   \begin{cases}
     \textrm{error} & \mbox{if } x\in \{\dots, -3, -2, -1, 0\}\\
     \Psi(x) & \mbox{if } x\not\in \{\dots, -3, -2, -1, 0\}\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
\f]
*
* @param x argument
* @return natural logarithm of the gamma function applied to
* argument
*/
inline double lgamma(double x) {
#if !__MINGW32__
  int sign = 1;
  return ::lgamma_r(x, &sign);
#else
  if (unlikely(x == 0.0))
    return std::numeric_limits<double>::infinity();
  return boost::math::lgamma(x, boost_policy_t());
#endif
}

/**
 * Return the natural logarithm of the gamma function applied
 * to the specified argument.
 *
 * @param x argument
 * @return natural logarithm of the gamma function applied to
 * argument
 */
inline double lgamma(int x) {
#if !__MINGW32__
  int sign = 1;
  return ::lgamma_r(x, &sign);
#else
  if (unlikely(x == 0.0))
    return std::numeric_limits<double>::infinity();
  return boost::math::lgamma(x, boost_policy_t());
#endif
}

}  // namespace math
}  // namespace stan
#endif
