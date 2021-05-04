#ifndef STAN_MATH_FWD_FUN_BETA_HPP
#define STAN_MATH_FWD_FUN_BETA_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/beta.hpp>
#include <stan/math/prim/fun/digamma.hpp>

namespace stan {
namespace math {

/**
 * Return fvar with the beta function applied to the specified
 * arguments and its gradient.
 *
 * The beta function is defined for \f$a > 0\f$ and \f$b > 0\f$ by
 *
 * \f$\mbox{B}(a, b) = \frac{\Gamma(a) \Gamma(b)}{\Gamma(a+b)}\f$.
 *
   \f[
   \mbox{beta}(\alpha, \beta) =
   \begin{cases}
     \int_0^1 u^{\alpha - 1} (1 - u)^{\beta - 1} \, du & \mbox{if } \alpha,
 \beta>0 \\[6pt] \textrm{NaN} & \mbox{if } \alpha = \textrm{NaN or } \beta =
 \textrm{NaN} \end{cases} \f]

   \f[
   \frac{\partial\, \mbox{beta}(\alpha, \beta)}{\partial \alpha} =
   \begin{cases}
     \left(\psi(\alpha)-\psi(\alpha+\beta)\right)*\mbox{beta}(\alpha, \beta)
      & \mbox{if } \alpha, \beta>0 \\[6pt]
     \textrm{NaN} & \mbox{if } \alpha = \textrm{NaN or } \beta = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{beta}(\alpha, \beta)}{\partial \beta} =
   \begin{cases}
     \left(\psi(\beta)-\psi(\alpha+\beta)\right)*\mbox{beta}(\alpha, \beta)
      & \mbox{if } \alpha, \beta>0 \\[6pt]
     \textrm{NaN} & \mbox{if } \alpha = \textrm{NaN or } \beta = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam T inner type of the fvar
 * @param x1 First value
 * @param x2 Second value
 * @return Fvar with result beta function of arguments and gradients.
 */
template <typename T>
inline fvar<T> beta(const fvar<T>& x1, const fvar<T>& x2) {
  const T beta_ab = beta(x1.val_, x2.val_);
  return fvar<T>(beta_ab,
                 beta_ab
                     * (x1.d_ * digamma(x1.val_) + x2.d_ * digamma(x2.val_)
                        - (x1.d_ + x2.d_) * digamma(x1.val_ + x2.val_)));
}

template <typename T>
inline fvar<T> beta(double x1, const fvar<T>& x2) {
  const T beta_ab = beta(x1, x2.val_);
  return fvar<T>(beta_ab,
                 x2.d_ * (digamma(x2.val_) - digamma(x1 + x2.val_)) * beta_ab);
}

template <typename T>
inline fvar<T> beta(const fvar<T>& x1, double x2) {
  const T beta_ab = beta(x1.val_, x2);
  return fvar<T>(beta_ab,
                 x1.d_ * (digamma(x1.val_) - digamma(x1.val_ + x2)) * beta_ab);
}

}  // namespace math
}  // namespace stan
#endif
