#ifndef STAN_MATH_REV_SCAL_FUN_ABS_HPP
#define STAN_MATH_REV_SCAL_FUN_ABS_HPP

#include <stan/math/rev/cplx.hpp>
#include <stan/math/rev/scal/fun/fabs.hpp>
#include <stan/math/rev/scal/fun/hypot.hpp>

namespace stan {
namespace math {

/**
 * Return the absolute value of the variable (std).
 *
 * Delegates to <code>fabs()</code> (see for doc).
 *
   \f[
   \mbox{abs}(x) =
   \begin{cases}
     |x| & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{abs}(x)}{\partial x} =
   \begin{cases}
     -1 & \mbox{if } x < 0 \\
     0 & \mbox{if } x = 0 \\
     1 & \mbox{if } x > 0 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Variable input.
 * @return Absolute value of variable.
 */
inline var abs(const var& a) { return fabs(a); }

}  // namespace math
}  // namespace stan

namespace std {
// constrained complex overload to forward zeroing var to std::abs
inline stan::math::var abs(std::complex<stan::math::var> const& t) {
  return abs(std::complex<stan::math::zeroing<stan::math::var>>(t));
}
}  // namespace std
#endif
