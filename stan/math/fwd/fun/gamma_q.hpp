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
 * The derivative with respect to the first argument is grad_reg_inc_gamma.
 * The derivative with respect to the second argument is minus the gamma
 * density, evaluated in log space so that neither pow nor tgamma can
 * overflow; this is the same form that gamma_p uses.
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
