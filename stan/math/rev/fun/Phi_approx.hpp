#ifndef STAN_MATH_REV_FUN_PHI_APPROX_HPP
#define STAN_MATH_REV_FUN_PHI_APPROX_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/inv_logit.hpp>
#include <stan/math/prim/fun/Phi_approx.hpp>

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
inline var Phi_approx(const var& a) {
  double av_squared = a.val() * a.val();
  double f = inv_logit(0.07056 * a.val() * av_squared + 1.5976 * a.val());
  double da = f * (1 - f) * (3.0 * 0.07056 * av_squared + 1.5976);
  return make_callback_var(
      f, [a, da](auto& vi) mutable { a.adj() += vi.adj() * da; });
}

template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto Phi_approx(const T& a) {
  arena_t<value_type_t<T>> f(a.rows(), a.cols());
  arena_t<value_type_t<T>> da(a.rows(), a.cols());
  for (Eigen::Index j = 0; j < a.cols(); ++j) {
    for (Eigen::Index i = 0; i < a.rows(); ++i) {
      const auto a_val = a.val().coeff(i, j);
      const auto av_squared = a_val * a_val;
      f.coeffRef(i, j) = inv_logit(0.07056 * a_val * av_squared
                                   + 1.5976 * a.val().coeff(i, j));
      da.coeffRef(i, j) = f.coeff(i, j) * (1 - f.coeff(i, j))
                          * (3.0 * 0.07056 * av_squared + 1.5976);
    }
  }
  return make_callback_var(f, [a, da](auto& vi) mutable {
    a.adj().array() += vi.adj().array() * da.array();
  });
}

}  // namespace math
}  // namespace stan
#endif
