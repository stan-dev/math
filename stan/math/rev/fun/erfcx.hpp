#ifndef STAN_MATH_REV_FUN_ERFCX_HPP
#define STAN_MATH_REV_FUN_ERFCX_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfcx.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * The scaled complementary error function for variables.
 *
 * The derivative is
 *
 * \f$\frac{d}{dx}\mbox{erfcx}(x) = 2x\,\mbox{erfcx}(x) -
 \frac{2}{\sqrt{\pi}}\f$
 *
 * which reuses the function value, so no extra `exp` or `erfc` evaluation is
 * needed and the derivative inherits the value's accuracy in the tails.
 *
   \f[
   \mbox{erfcx}(x) =
   \begin{cases}
     \exp(x^2)\operatorname{erfc}(x) & \mbox{if } -\infty\leq x \leq \infty
     \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{erfcx}(x)}{\partial x} =
   \begin{cases}
     2x\operatorname{erfcx}(x) - \frac{2}{\sqrt{\pi}} & \mbox{if }
     -\infty\leq x\leq \infty \\[6pt] \textrm{NaN} & \mbox{if }
     x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a The variable.
 * @return Scaled complementary error function applied to the variable.
 */
inline var erfcx(const var& a) {
  double val = erfcx(a.val());
  return make_callback_var(val, [a, val](auto& vi) mutable {
    a.adj() += vi.adj() * (2.0 * a.val() * val - TWO_OVER_SQRT_PI);
  });
}

/**
 * The scaled complementary error function for matrix variables.
 *
 * @tparam T a matrix type
 * @param a The variable.
 * @return Scaled complementary error function applied elementwise.
 */
template <typename T, require_matrix_t<T>* = nullptr>
inline auto erfcx(const var_value<T>& a) {
  auto val = to_arena(erfcx(a.val()));
  return make_callback_var(val, [a, val](auto& vi) mutable {
    a.adj().array()
        += vi.adj().array()
           * (2.0 * a.val().array() * val.array() - TWO_OVER_SQRT_PI);
  });
}

}  // namespace math
}  // namespace stan
#endif
