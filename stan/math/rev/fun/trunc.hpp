#ifndef STAN_MATH_REV_FUN_TRUNC_HPP
#define STAN_MATH_REV_FUN_TRUNC_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_nan.hpp>
#include <stan/math/prim/fun/trunc.hpp>

namespace stan {
namespace math {

/**
 * Returns the truncation of the specified variable (C99).
 *
 * See ::trunc() for the double-based version.
 *
 * The derivative is zero everywhere but at integer values, so for
 * convenience the derivative is defined to be everywhere zero,
 *
 * \f$\frac{d}{dx} \mbox{trunc}(x) = 0\f$.
 *
 *
   \f[
   \mbox{trunc}(x) =
   \begin{cases}
     \lfloor x \rfloor & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{trunc}(x)}{\partial x} =
   \begin{cases}
     0 & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Specified variable.
 * @return Truncation of the variable.
 */
inline var trunc(const var& a) { return var(trunc(a.val())); }

}  // namespace math
}  // namespace stan
#endif
