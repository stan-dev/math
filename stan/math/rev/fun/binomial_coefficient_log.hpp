#ifndef STAN_MATH_REV_FUN_BINOMIAL_COEFFICIENT_LOG_HPP
#define STAN_MATH_REV_FUN_BINOMIAL_COEFFICIENT_LOG_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/binomial_coefficient_log.hpp>
#include <stan/math/prim/fun/digamma.hpp>

namespace stan {
namespace math {

namespace internal {
class binomial_coefficient_log_vv_vari : public op_vv_vari {
 public:
  binomial_coefficient_log_vv_vari(vari* avi, vari* bvi)
      : op_vv_vari(binomial_coefficient_log(avi->val_, bvi->val_), avi, bvi) {}
  void chain() {
    double digamma_ambp1 = digamma(avi_->val_ - bvi_->val_ + 1);

    avi_->adj_ += adj_ * (digamma(avi_->val_ + 1) - digamma_ambp1);
    bvi_->adj_ += adj_ * (digamma_ambp1 - digamma(bvi_->val_ + 1));
  }
};

class binomial_coefficient_log_vd_vari : public op_vd_vari {
 public:
  binomial_coefficient_log_vd_vari(vari* avi, double b)
      : op_vd_vari(binomial_coefficient_log(avi->val_, b), avi, b) {}
  void chain() {
    avi_->adj_ += adj_ * (digamma(avi_->val_ + 1) - digamma(avi_->val_ - bd_ + 1));
  }
};

class binomial_coefficient_log_dv_vari : public op_dv_vari {
 public:
  binomial_coefficient_log_dv_vari(double a, vari* bvi)
      : op_dv_vari(binomial_coefficient_log(a, bvi->val_), a, bvi) {}
  void chain() {
    bvi_->adj_ += adj_ * (digamma(ad_ - bvi_->val_ + 1) - digamma(bvi_->val_ + 1));
  }
};
}  // namespace internal

/**
 * Return the log of the binomial coefficient for the specified
 * arguments and its gradients.
 * 
 * See the docs for the prim version for all relevant formulae.
 * @param a var Argument
 * @param b var Argument
 * @return Result of log (a choose b)
 */
inline var binomial_coefficient_log(const var& a, const var& b) {
  return var(new internal::binomial_coefficient_log_vv_vari(a.vi_, b.vi_));
}

/**
 * Return the log of the binomial coefficient for the specified
 * arguments and its gradients.
 * 
 * See the docs for the prim version for all relevant formulae.
 * @param a var Argument
 * @param b double Argument
 * @return Result of log (a choose b)
 */
inline var binomial_coefficient_log(const var& a, double b) {
  return var(new internal::binomial_coefficient_log_vd_vari(a.vi_, b));
}

/**
 * Return the log of the binomial coefficient for the specified
 * arguments and its gradients.
 * 
 * See the docs for the prim version for all relevant formulae.
 * @param a double Argument
 * @param b var Argument
 * @return Result of log (a choose b)
 */
inline var binomial_coefficient_log(double a, const var& b) {
  return var(new internal::binomial_coefficient_log_dv_vari(a, b.vi_));
}

}  // namespace math
}  // namespace stan
#endif
