#ifndef STAN_MATH_REV_FUN_INV_SQRT_HPP
#define STAN_MATH_REV_FUN_INV_SQRT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/inv_sqrt.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 *
   \f[
   \mbox{inv\_sqrt}(x) =
   \begin{cases}
     \frac{1}{\sqrt{x}} & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{inv\_sqrt}(x)}{\partial x} =
   \begin{cases}
     -\frac{1}{2\sqrt{x^3}} & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 */
inline var inv_sqrt(const var& a) {
  auto denom = a.val() * std::sqrt(a.val());
  return make_callback_var(inv_sqrt(a.val()), [a, denom](auto& vi) mutable {
    a.adj() -= 0.5 * vi.adj() / denom;
  });
}

}  // namespace math
}  // namespace stan
#endif
