#ifndef STAN_MATH_FWD_FUN_LOG_INV_LOGIT_DIFF_HPP
#define STAN_MATH_FWD_FUN_LOG_INV_LOGIT_DIFF_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/inv_logit.hpp>
#include <stan/math/fwd/fun/inv.hpp>
#include <stan/math/fwd/fun/expm1.hpp>
#include <stan/math/prim/fun/log_inv_logit_diff.hpp>

namespace stan {
namespace math {

/**
 * Returns fvar with the natural logarithm of the difference of the
 * inverse logits of the specified arguments and its gradients.
 *
   \f[
     \mathrm{log\_inv\_logit\_diff}(x,y) =
      \ln\left(\frac{1}{1+\exp(-x)}-\frac{1}{1+\exp(-y)}\right)
   \f]

   \f[
    \frac{\partial}{\partial x} = \mbox{logit}^{-1}(-x) + \frac{1}{e^{x-y}-1}
   \f]

   \f[
    \frac{\partial}{\partial y} = -\mbox{logit}^{-1}(y) - \frac{1}{e^{x-y}-1}
   \f]
 *
 * @tparam T inner type of the fvar
 * @param x Argument.
 * @param y Argument.
 * @return Fvar with result of log difference of inverse logits of arguments
 *           and gradients.
 */
template <typename T>
inline fvar<T> log_inv_logit_diff(const fvar<T>& x, const fvar<T>& y) {
  // both partials as sums of positive terms, no cancellation for x >> 0
  const auto c = inv(expm1(x.val_ - y.val_));
  return fvar<T>(
      log_inv_logit_diff(x.val_, y.val_),
      x.d_ * (c + inv_logit(-x.val_)) - y.d_ * (c + inv_logit(y.val_)));
}

template <typename T>
inline fvar<T> log_inv_logit_diff(const fvar<T>& x, double y) {
  return fvar<T>(log_inv_logit_diff(x.val_, y),
                 x.d_ * (inv(expm1(x.val_ - y)) + inv_logit(-x.val_)));
}

template <typename T>
inline fvar<T> log_inv_logit_diff(double x, const fvar<T>& y) {
  return fvar<T>(log_inv_logit_diff(x, y.val_),
                 -y.d_ * (inv(expm1(x - y.val_)) + inv_logit(y.val_)));
}

}  // namespace math
}  // namespace stan
#endif
