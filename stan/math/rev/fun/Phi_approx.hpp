#ifndef STAN_MATH_REV_FUN_PHI_APPROX_HPP
#define STAN_MATH_REV_FUN_PHI_APPROX_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>

namespace stan {
namespace math {

/**
 * Approximation of the unit normal CDF for variables (stan).
 *
 * http://www.jiem.org/index.php/jiem/article/download/60/27
 *
 *
   \f[
   \mbox{Phi\_approx}(x) =
   \begin{cases}
     \Phi_{\mbox{\footnotesize approx}}(x) & \mbox{if } -\infty\leq x\leq \infty
 \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN} \end{cases} \f]

   \f[
   \frac{\partial\, \mbox{Phi\_approx}(x)}{\partial x} =
   \begin{cases}
     \frac{\partial\, \Phi_{\mbox{\footnotesize approx}}(x)}{\partial x}
     & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \Phi_{\mbox{\footnotesize approx}}(x) = \mbox{logit}^{-1}(0.07056 \,
   x^3 + 1.5976 \, x)
   \f]

   \f[
   \frac{\partial \, \Phi_{\mbox{\footnotesize approx}}(x)}{\partial x}
   = -\Phi_{\mbox{\footnotesize approx}}^2(x)
   e^{-0.07056x^3 - 1.5976x}(-0.21168x^2-1.5976)
   \f]
 *
 * @param a Variable argument.
 * @return The corresponding unit normal cdf approximation.
 */
template <typename T>
inline var_value<T> Phi_approx(const var_value<T>& a) {
  auto av = a.vi_->val_;
  auto av_squared = av * av;
  auto av_cubed = av * av_squared;
  auto f = inv_logit(0.07056 * av_cubed + 1.5976 * av);
  auto da = f * (1 - f) * (3.0 * 0.07056 * av_squared + 1.5976);
  return {new precomp_v_vari<T, vari_value<T>>(f, a.vi_, da)};
}

}  // namespace math
}  // namespace stan
#endif
