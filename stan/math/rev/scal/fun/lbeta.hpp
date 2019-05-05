#ifndef STAN_MATH_REV_SCAL_FUN_LBETA_HPP
#define STAN_MATH_REV_SCAL_FUN_LBETA_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/scal/fun/lbeta.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>

namespace stan {
namespace math {

/*
 * Returns the natural logarithm of the beta function and its gradients.
 *
   \f[
     \mathrm{lbeta}(a,b) = \ln\left(B\left(a,b\right)\right)
   \f]

   \f[
    \frac{\partial }{\partial a} = \psi^{\left(0\right)}\left(a\right)
                                      - \psi^{\left(0\right)}\left(a + b\right)
   \f]

   \f[
    \frac{\partial }{\partial b} = \psi^{\left(0\right)}\left(b\right)
                                      - \psi^{\left(0\right)}\left(a + b\right)
   \f]
 *
 * @tparam T1 Type of a argument
 * @tparam T2 Type of b argument
 * @param a Argument
 * @param b Argument
 * @return Result of log beta function
 */
namespace internal {
class lbeta_vv_vari : public op_vv_vari {
 public:
  lbeta_vv_vari(vari* avi, vari* bvi)
      : op_vv_vari(lbeta(avi->val_, bvi->val_), avi, bvi) {}
  void chain() {
    double digamma_ab = digamma(avi_->val_ + bvi_->val_);
    avi_->adj_ += adj_ * digamma(avi_->val_) - digamma_ab;

    bvi_->adj_ += adj_ * digamma(bvi_->val_) - digamma_ab;
  }
};

class lbeta_vd_vari : public op_vd_vari {
 public:
  lbeta_vd_vari(vari* avi, double b)
      : op_vd_vari(lbeta(avi->val_, b), avi, b) {}
  void chain() {
    avi_->adj_ += adj_ * digamma(avi_->val_) - digamma(avi_->val_ + bd_);
  }
};

class lbeta_dv_vari : public op_dv_vari {
 public:
  lbeta_dv_vari(double a, vari* bvi)
      : op_dv_vari(lbeta(a, bvi->val_), a, bvi) {}
  void chain() {
    bvi_->adj_ += adj_ * digamma(bvi_->val_) - digamma(ad_ + bvi_->val_);
  }
};
}  // namespace internal

inline var lbeta(const var& a, double b) {
  return var(new internal::lbeta_vd_vari(a.vi_, b));
}

inline var lbeta(const var& a, const var& b) {
  return var(new internal::lbeta_vv_vari(a.vi_, b.vi_));
}

inline var lbeta(double a, const var& b) {
  return var(new internal::lbeta_dv_vari(a, b.vi_));
}

}  // namespace math
}  // namespace stan
#endif
