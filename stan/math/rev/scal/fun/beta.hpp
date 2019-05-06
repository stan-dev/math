#ifndef STAN_MATH_REV_SCAL_FUN_BETA_HPP
#define STAN_MATH_REV_SCAL_FUN_BETA_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/scal/fun/beta.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>

namespace stan {
namespace math {

/*
 * Returns the beta function and its gradients.
 *
   \f[
     \mathrm{beta}(a,b) = \left(B\left(a,b\right)\right)
   \f]

   \f[
    \frac{\partial }{\partial a} = \left(\psi^{\left(0\right)}\left(a\right)
                                      - \psi^{\left(0\right)}\left(a + b\right)\right)
                                    * \mathrm{beta}(a,b)
   \f]

   \f[
    \frac{\partial }{\partial b} = \left(\psi^{\left(0\right)}\left(b\right)
                                      - \psi^{\left(0\right)}\left(a + b\right)\right)
                                    * \mathrm{beta}(a,b)
   \f]
 *
 * @tparam T1 Type of a argument
 * @tparam T2 Type of b argument
 * @param a Argument
 * @param b Argument
 * @return Result of beta function
 */
namespace internal {
class beta_vv_vari : public op_vv_vari {
 public:
  beta_vv_vari(vari* avi, vari* bvi)
      : op_vv_vari(beta(avi->val_, bvi->val_), avi, bvi) {}
  void chain() {
    const double beta_ab = beta(avi_->val_, bvi_->val_);
    const double digamma_ab = digamma(avi_->val_ + bvi_->val_);
    avi_->adj_ += adj_ * (digamma(avi_->val_) - digamma_ab) * beta_ab;

    bvi_->adj_ += adj_ * (digamma(bvi_->val_) - digamma_ab) * beta_ab;
  }
};

class beta_vd_vari : public op_vd_vari {
 public:
  beta_vd_vari(vari* avi, double b)
      : op_vd_vari(beta(avi->val_, b), avi, b) {}
  void chain() {
    avi_->adj_ += adj_ * (digamma(avi_->val_) - digamma(avi_->val_ + bd_))
                          * beta(avi_->val_, bd_);
  }
};

class beta_dv_vari : public op_dv_vari {
 public:
  beta_dv_vari(double a, vari* bvi)
      : op_dv_vari(beta(a, bvi->val_), a, bvi) {}
  void chain() {
    bvi_->adj_ += adj_ * (digamma(bvi_->val_) - digamma(ad_ + bvi_->val_))
                          * beta(ad_, bvi_->val_);
  }
};
}  // namespace internal

inline var beta(const var& a, double b) {
  return var(new internal::beta_vd_vari(a.vi_, b));
}

inline var beta(const var& a, const var& b) {
  return var(new internal::beta_vv_vari(a.vi_, b.vi_));
}

inline var beta(double a, const var& b) {
  return var(new internal::beta_dv_vari(a, b.vi_));
}

}  // namespace math
}  // namespace stan
#endif
