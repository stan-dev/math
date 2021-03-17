#ifndef STAN_MATH_REV_FUN_INV_SQUARE_HPP
#define STAN_MATH_REV_FUN_INV_SQUARE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/inv_square.hpp>

namespace stan {
namespace math {

/**
 *
   \f[
   \mbox{inv\_square}(x) =
   \begin{cases}
     \frac{1}{x^2} & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{inv\_square}(x)}{\partial x} =
   \begin{cases}
     -\frac{2}{x^3} & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 */
inline var inv_square(const var& a) {
  auto a_cube = a.val() * a.val() * a.val();
  return make_callback_var(inv_square(a.val()), [a, a_cube](auto& vi) mutable {
    a.adj() -= 2 * vi.adj() / a_cube;
  });
}

}  // namespace math
}  // namespace stan
#endif
