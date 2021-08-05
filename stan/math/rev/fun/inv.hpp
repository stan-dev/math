#ifndef STAN_MATH_REV_FUN_INV_HPP
#define STAN_MATH_REV_FUN_INV_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/inv.hpp>

namespace stan {
namespace math {

/**
 *
   \f[
   \mbox{inv}(x) =
   \begin{cases}
     \frac{1}{x} & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{inv}(x)}{\partial x} =
   \begin{cases}
     -\frac{1}{x^2} & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 */
inline var inv(const var& a) {
  return make_callback_var(inv(a.val()), [a](auto& vi) mutable {
    a.adj() -= vi.adj() / (a.val() * a.val());
  });
}

}  // namespace math
}  // namespace stan
#endif
