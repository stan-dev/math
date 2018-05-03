#ifndef STAN_MATH_FWD_SCAL_FUN_LOG_INV_LOGIT_DIFF_HPP
#define STAN_MATH_FWD_SCAL_FUN_LOG_INV_LOGIT_DIFF_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/scal/fun/log_inv_logit_diff.hpp>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> log_inv_logit_diff(const fvar<T>& x, const fvar<T>& y) {
  return fvar<T>(
      log_inv_logit_diff(x.val_, y.val_),
      x.d_ * (-exp(x.val_)/(exp(y.val_) - exp(x.val_))
               - exp(x.val_)/(exp(x.val_) + 1.0))
        + y.d_ * (-exp(y.val_)/(exp(x.val_) - exp(y.val_))
                          - exp(y.val_)/(exp(y.val_) + 1.0)));
}

template <typename T>
inline fvar<T> log_inv_logit_diff(const fvar<T>& x, double y) {
  return fvar<T>(log_inv_logit_diff(x.val_, y),
                 x.d_ * (-exp(x.val_)/(exp(y) - exp(x.val_))
                          - exp(x.val_)/(exp(x.val_) + 1.0)));
}

template <typename T>
inline fvar<T> log_inv_logit_diff(double x, const fvar<T>& y) {
  return fvar<T>(log_inv_logit_diff(x, y.val_),
                 y.d_ * (-exp(y.val_)/(exp(x) - exp(y.val_))
                          - exp(y.val_)/(exp(y.val_) + 1.0)));
}

}  // namespace math
}  // namespace stan
#endif
