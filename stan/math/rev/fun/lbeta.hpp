#ifndef STAN_MATH_REV_FUN_LBETA_HPP
#define STAN_MATH_REV_FUN_LBETA_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/digamma.hpp>

namespace stan {
namespace math {

namespace internal {
class lbeta_vv_vari : public op_vv_vari {
 public:
  lbeta_vv_vari(vari* avi, vari* bvi)
      : op_vv_vari(lbeta(avi->val_, bvi->val_), avi, bvi) {}
  void chain() {
    const double digamma_ab = digamma(avi_->val_ + bvi_->val_);
    avi_->adj_ += adj_ * (digamma(avi_->val_) - digamma_ab);

    bvi_->adj_ += adj_ * (digamma(bvi_->val_) - digamma_ab);
  }
};

class lbeta_vd_vari : public op_vd_vari {
 public:
  lbeta_vd_vari(vari* avi, double b)
      : op_vd_vari(lbeta(avi->val_, b), avi, b) {}
  void chain() {
    avi_->adj_ += adj_ * (digamma(avi_->val_) - digamma(avi_->val_ + bd_));
  }
};

class lbeta_dv_vari : public op_dv_vari {
 public:
  lbeta_dv_vari(double a, vari* bvi)
      : op_dv_vari(lbeta(a, bvi->val_), a, bvi) {}
  void chain() {
    bvi_->adj_ += adj_ * (digamma(bvi_->val_) - digamma(ad_ + bvi_->val_));
  }
};
}  // namespace internal

/**
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
 * @param a var Argument
 * @param b var Argument
 * @return Result of log beta function
 */
inline var lbeta(const var& a, const var& b) {
  return var(new internal::lbeta_vv_vari(a.vi_, b.vi_));
}

/**
 * Returns the natural logarithm of the beta function and its gradients.
 *
   \f[
     \mathrm{lbeta}(a,b) = \ln\left(B\left(a,b\right)\right)
   \f]

   \f[
    \frac{\partial }{\partial a} = \psi^{\left(0\right)}\left(a\right)
                                      - \psi^{\left(0\right)}\left(a + b\right)
   \f]
 * @param a var Argument
 * @param b double Argument
 * @return Result of log beta function
 */
inline var lbeta(const var& a, double b) {
  return var(new internal::lbeta_vd_vari(a.vi_, b));
}

/**
 * Returns the natural logarithm of the beta function and its gradients.
 *
   \f[
     \mathrm{lbeta}(a,b) = \ln\left(B\left(a,b\right)\right)
   \f]

   \f[
    \frac{\partial }{\partial b} = \psi^{\left(0\right)}\left(b\right)
                                      - \psi^{\left(0\right)}\left(a + b\right)
   \f]
 * @param a double Argument
 * @param b var Argument
 * @return Result of log beta function
 */
inline var lbeta(double a, const var& b) {
  return var(new internal::lbeta_dv_vari(a, b.vi_));
}

}  // namespace math
}  // namespace stan
#endif
