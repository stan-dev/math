#ifndef STAN_MATH_REV_FUN_ERF_HPP
#define STAN_MATH_REV_FUN_ERF_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * The error function for variables (C99).
 *
 * The derivative is
 *
 * \f$\frac{d}{dx} \mbox{erf}(x) = \frac{2}{\sqrt{\pi}} \exp(-x^2)\f$.
 *
 *
   \f[
   \mbox{erf}(x) =
   \begin{cases}
     \operatorname{erf}(x) & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{erf}(x)}{\partial x} =
   \begin{cases}
     \frac{\partial\, \operatorname{erf}(x)}{\partial x} & \mbox{if }
 -\infty\leq x\leq \infty \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \operatorname{erf}(x)=\frac{2}{\sqrt{\pi}}\int_0^x e^{-t^2}dt
   \f]

   \f[
   \frac{\partial \, \operatorname{erf}(x)}{\partial x} = \frac{2}{\sqrt{\pi}}
 e^{-x^2} \f]
 *
 * @param a The variable.
 * @return Error function applied to the variable.
 */
inline var erf(const var& a) {
  auto precomp_erf = TWO_OVER_SQRT_PI * std::exp(-a.val() * a.val());
  return make_callback_var(erf(a.val()), [a, precomp_erf](auto& vi) mutable {
    a.adj() += vi.adj() * precomp_erf;
  });
}

template <typename T, require_matrix_t<T>* = nullptr>
inline auto erf(const var_value<T>& a) {
  auto precomp_erf
      = to_arena(TWO_OVER_SQRT_PI * (-a.val().array().square()).exp());
  return make_callback_var(erf(a.val()), [a, precomp_erf](auto& vi) mutable {
    a.adj().array() += vi.adj().array() * precomp_erf;
  });
}

}  // namespace math
}  // namespace stan
#endif
