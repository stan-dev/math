#ifndef STAN_MATH_REV_FUN_LOG_ADD_EXP_HPP
#define STAN_MATH_REV_FUN_LOG_ADD_EXP_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/rev/fun/inv_logit.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log_add_exp.hpp>
#include <cmath>
#include <vector>

#include <iostream>

namespace stan {
namespace math {
namespace internal {

class log_add_exp_vv_vari : public op_vv_vari {
 public:
  log_add_exp_vv_vari(vari* avi, vari* bvi)
      : op_vv_vari(log_add_exp(avi->val_, bvi->val_), avi, bvi) {}
  void chain() {
    double exp_a = std::exp(avi_->val_);
    double exp_b = std::exp(bvi_->val_);
    double sum_exp = exp_a + exp_b;

    avi_->adj_ += adj_ * (exp_a / sum_exp);
    bvi_->adj_ += adj_ * (exp_b / sum_exp);
  }
};

class log_add_exp_vd_vari : public op_vd_vari {
 public:
  log_add_exp_vd_vari(vari* avi, double b)
      : op_vd_vari(log_add_exp(avi->val_, b), avi, b) {}
  void chain() {
    if (val_ == NEGATIVE_INFTY) {
      avi_->adj_ += adj_;
    } else {
      double exp_a = std::exp(avi_->val_);
      avi_->adj_ += adj_ * (exp_a / (exp_a + std::exp(bd_)));
    }
  }
};

}  // namespace internal

/**
 * Returns the element-wise log sum of exponentials.
 */
inline var log_add_exp(const var& a, const var& b) {
  return var(new internal::log_add_exp_vv_vari(a.vi_, b.vi_));
}

/**
 * Returns the log sum of exponentials.
 */
inline var log_add_exp(const var& a, double b) {
  return var(new internal::log_add_exp_vd_vari(a.vi_, b));
}

/**
 * Returns the element-wise log sum of exponentials.
 */
inline var log_add_exp(double a, const var& b) {
  return var(new internal::log_add_exp_vd_vari(b.vi_, a));
}

/**
 * Returns element-wise log sum of exponentials for Eigen types.
 *
 * @tparam T A type inheriting from EigenBase with var scalar type
 * @param x First input
 * @param y Second input
 */
template <typename T, require_eigen_st<is_var, T>* = nullptr>
inline T log_add_exp(const T& x, const T& y) {
  return apply_scalar_binary(
      x, y, [](const auto& a, const auto& b) { return log_add_exp(a, b); });
}

}  // namespace math
}  // namespace stan

#endif
