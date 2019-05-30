#ifndef STAN_MATH_PRIM_SCAL_FUN_LGAMMA_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LGAMMA_HPP

#include <limits>
#include <boost/math/special_functions/gamma.hpp>
#include <stan/math/prim/scal/meta/likely.hpp>
#include <stan/math/prim/scal/fun/boost_policy.hpp>

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
* Note: The use of std::lgamma is explicitly avoided here. The std
* implementation is not suitable for concurrent use in threaded
* programs as the exact behavior under threading is implementation
* specific. See discussion under https://github.com/stan-dev/math/issues/1250
*
* @param x argument
* @return natural logarithm of the gamma function applied to
* argument
*/
inline double lgamma(double x) {
  if (unlikely(x == 0.0))
    return std::numeric_limits<double>::infinity();
  return boost::math::lgamma(x, boost_policy_t());
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
  if (unlikely(x == 0))
    return std::numeric_limits<double>::infinity();
  return boost::math::lgamma(x, boost_policy_t());
}

}  // namespace math
}  // namespace stan
#endif
