#ifndef STAN_MATH_FWD_FUN_GAMMA_Q_HPP
#define STAN_MATH_FWD_FUN_GAMMA_Q_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/digamma.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <stan/math/fwd/fun/lgamma.hpp>
#include <stan/math/fwd/fun/log.hpp>
#include <stan/math/fwd/fun/tgamma.hpp>
#include <stan/math/prim/fun/gamma_q.hpp>
#include <stan/math/prim/fun/grad_reg_inc_gamma.hpp>
#include <cmath>

namespace stan {
namespace math {

/*
 * Two changes from the previous version of this file.
 *
 * The derivative with respect to the shape parameter used to be an inlined
 * copy of the series in `grad_reg_inc_gamma`, with the same hard-coded 1e-6
 * tolerance and without that function's second branch. The copy is now
 * replaced by a call to the root itself.
 *
 * Keeping the copy hid the defect it shared with the root. In
 * `test/prob/chi_square`, the generated `ffv` case compares the
 * distribution's analytic partials against autodiff through
 * `log(gamma_q(nu * 0.5, y * 0.5))`. Both routes used the same inaccurate
 * series, so their errors cancelled and the comparison passed. Once the
 * root became accurate and the copy did not, the same comparison failed by
 * 2.4e-03 at third order.
 *
 * The derivative with respect to the second argument used to be
 * `-exp(-x2) * pow(x2, x1 - 1) / tgamma(x1)`. Both `pow` and `tgamma`
 * overflow for moderate arguments, and the quotient then silently returns
 * -0 rather than a small negative number. It is now evaluated in log space,
 * which is the form `gamma_p` already uses for the same quantity.
 */

template <typename T>
inline fvar<T> gamma_q(const fvar<T>& x1, const fvar<T>& x2) {
  T u = gamma_q(x1.val_, x2.val_);

  T der1
      = grad_reg_inc_gamma(x1.val_, x2.val_, tgamma(x1.val_), digamma(x1.val_));
  T der2 = -exp(-x2.val_ + (x1.val_ - 1.0) * log(x2.val_) - lgamma(x1.val_));

  return fvar<T>(u, x1.d_ * der1 + x2.d_ * der2);
}

template <typename T>
inline fvar<T> gamma_q(const fvar<T>& x1, double x2) {
  T u = gamma_q(x1.val_, x2);

  T der1 = grad_reg_inc_gamma(x1.val_, x2, tgamma(x1.val_), digamma(x1.val_));

  return fvar<T>(u, x1.d_ * der1);
}

template <typename T>
inline fvar<T> gamma_q(double x1, const fvar<T>& x2) {
  T u = gamma_q(x1, x2.val_);

  T der2 = -exp(-x2.val_ + (x1 - 1.0) * log(x2.val_) - lgamma(x1));

  return fvar<T>(u, x2.d_ * der2);
}
}  // namespace math
}  // namespace stan
#endif
