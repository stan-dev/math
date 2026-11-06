#ifndef STAN_MATH_REV_FUN_LBETA_HPP
#define STAN_MATH_REV_FUN_LBETA_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/digamma.hpp>
#include <stan/math/prim/fun/digamma_diff.hpp>
#include <stan/math/prim/fun/lbeta.hpp>

namespace stan {
namespace math {

namespace internal {
/**
 * Return digamma(x) - digamma(x + y), the partial derivative of
 * lbeta(x, y) in x. For large x the plain difference loses all digits,
 * so from x = digamma_diff_min_x on it is formed by
 * digamma_diff. Below that value (and for NaN) the plain difference is
 * accurate enough for a gradient and cheaper.
 */
inline double lbeta_partial(double x, double y) {
  if (x >= digamma_diff_min_x) {
    return -digamma_diff(x, y);
  }
  return digamma(x) - digamma(x + y);
}

class lbeta_vv_vari : public op_vv_vari {
 public:
  lbeta_vv_vari(vari* avi, vari* bvi)
      : op_vv_vari(lbeta(avi->val_, bvi->val_), avi, bvi) {}
  void chain() {
    const double a = avi_->val_;
    const double b = bvi_->val_;
    if (a >= digamma_diff_min_x || b >= digamma_diff_min_x) {
      avi_->adj_ += adj_ * lbeta_partial(a, b);
      bvi_->adj_ += adj_ * lbeta_partial(b, a);
    } else {
      // both arguments small: the plain differences share digamma(a + b)
      const double digamma_ab = digamma(a + b);
      avi_->adj_ += adj_ * (digamma(a) - digamma_ab);
      bvi_->adj_ += adj_ * (digamma(b) - digamma_ab);
    }
  }
};

class lbeta_vd_vari : public op_vd_vari {
 public:
  lbeta_vd_vari(vari* avi, double b)
      : op_vd_vari(lbeta(avi->val_, b), avi, b) {}
  void chain() { avi_->adj_ += adj_ * lbeta_partial(avi_->val_, bd_); }
};

class lbeta_dv_vari : public op_dv_vari {
 public:
  lbeta_dv_vari(double a, vari* bvi)
      : op_dv_vari(lbeta(a, bvi->val_), a, bvi) {}
  void chain() { bvi_->adj_ += adj_ * lbeta_partial(bvi_->val_, ad_); }
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
