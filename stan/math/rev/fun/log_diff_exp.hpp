#ifndef STAN_MATH_REV_FUN_LOG_DIFF_EXP_HPP
#define STAN_MATH_REV_FUN_LOG_DIFF_EXP_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/expm1.hpp>
#include <stan/math/prim/fun/log_diff_exp.hpp>

namespace stan {
namespace math {

namespace internal {
class log_diff_exp_vv_vari : public op_vv_vari {
 public:
  log_diff_exp_vv_vari(vari* avi, vari* bvi)
      : op_vv_vari(log_diff_exp(avi->val_, bvi->val_), avi, bvi) {}
  void chain() {
    avi_->adj_ -= adj_ / expm1(bvi_->val_ - avi_->val_);
    bvi_->adj_ -= adj_ / expm1(avi_->val_ - bvi_->val_);
  }
};
class log_diff_exp_vd_vari : public op_vd_vari {
 public:
  log_diff_exp_vd_vari(vari* avi, double b)
      : op_vd_vari(log_diff_exp(avi->val_, b), avi, b) {}
  void chain() {
    if (val_ == NEGATIVE_INFTY) {
      avi_->adj_ += (bd_ == NEGATIVE_INFTY) ? adj_ : adj_ * INFTY;
    } else {
      avi_->adj_ -= adj_ / expm1(bd_ - avi_->val_);
    }
  }
};
class log_diff_exp_dv_vari : public op_dv_vari {
 public:
  log_diff_exp_dv_vari(double a, vari* bvi)
      : op_dv_vari(log_diff_exp(a, bvi->val_), a, bvi) {}
  void chain() {
    if (val_ == NEGATIVE_INFTY) {
      bvi_->adj_ -= adj_ * INFTY;
    } else {
      bvi_->adj_ -= adj_ / expm1(ad_ - bvi_->val_);
    }
  }
};
}  // namespace internal

/**
 * Returns the log difference of the exponentiated arguments.
 *
 * @param[in] a First argument.
 * @param[in] b Second argument.
 * @return Log difference of the exponentiated arguments.
 */
inline var log_diff_exp(const var& a, const var& b) {
  return var(new internal::log_diff_exp_vv_vari(a.vi_, b.vi_));
}

/**
 * Returns the log difference of the exponentiated arguments.
 *
 * @param[in] a First argument.
 * @param[in] b Second argument.
 * @return Log difference of the exponentiated arguments.
 */
inline var log_diff_exp(const var& a, double b) {
  return var(new internal::log_diff_exp_vd_vari(a.vi_, b));
}

/**
 * Returns the log difference of the exponentiated arguments.
 *
 * @param[in] a First argument.
 * @param[in] b Second argument.
 * @return Log difference of the exponentiated arguments.
 */
inline var log_diff_exp(double a, const var& b) {
  return var(new internal::log_diff_exp_dv_vari(a, b.vi_));
}

}  // namespace math
}  // namespace stan
#endif
