#ifndef STAN_MATH_REV_SCAL_FUN_GAMMA_P_HPP
#define STAN_MATH_REV_SCAL_FUN_GAMMA_P_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>
#include <stan/math/prim/scal/fun/gamma_p.hpp>
#include <stan/math/prim/scal/fun/tgamma.hpp>
#include <stan/math/prim/scal/fun/grad_reg_lower_inc_gamma.hpp>
#include <valarray>
#include <iostream>

namespace stan {
  namespace math {

    namespace {
      class gamma_p_vv_vari : public op_vv_vari {
      public:
        gamma_p_vv_vari(vari* avi, vari* bvi) :
          op_vv_vari(gamma_p(avi->val_, bvi->val_),
                     avi, bvi) {
        }
        void chain() {
          using std::fabs;
          using std::log;
          using std::exp;
          using std::pow;
          using boost::math::lgamma;
          // return zero derivative as gamma_p is flat
          // to machine precision for b / a > 10
          if (std::fabs(bvi_->val_ / avi_->val_) > 10 ) return;

          avi_->adj_ += adj_
            * grad_reg_lower_inc_gamma(avi_->val_, bvi_->val_, 1.0e-10);
          bvi_->adj_ += adj_ * exp(
            - bvi_->val_ + (avi_->val_ - 1.0) * log(bvi_->val_) 
            - lgamma(avi_->val_));

          std::cout << "Called it." << std::endl;
          std::cout << "adj_: " << adj_ << std::endl;
          std::cout << "avi_->val_: " << avi_->val_ << std::endl;
          std::cout << "bvi_->val_: " << bvi_->val_ << std::endl;
          std::cout << "lgamma(avi_->val_): " << lgamma(avi_->val_) << std::endl;

        }
      };

      class gamma_p_vd_vari : public op_vd_vari {
      public:
        gamma_p_vd_vari(vari* avi, double b) :
          op_vd_vari(gamma_p(avi->val_, b),
                     avi, b) {
        }
        void chain() {
          // return zero derivative as gamma_p is flat
          // to machine precision for b / a > 10
          if (std::fabs(bd_ / avi_->val_) > 10)
            return;

          avi_->adj_ += adj_
            * grad_reg_lower_inc_gamma(avi_->val_, bd_, 1.0e-10);
        }
      };

      class gamma_p_dv_vari : public op_dv_vari {
      public:
        gamma_p_dv_vari(double a, vari* bvi) :
          op_dv_vari(gamma_p(a, bvi->val_),
                     a, bvi) {
        }
        void chain() {
          // return zero derivative as gamma_p is flat to
          // machine precision for b / a > 10
          if (std::fabs(bvi_->val_ / ad_) > 10 )
            return;
          bvi_->adj_ += adj_
            * (std::exp(-bvi_->val_) * std::pow(bvi_->val_, ad_ - 1.0)
               / tgamma(ad_));
        }
      };
    }

    inline var gamma_p(const var& a,
                       const var& b) {
      return var(new gamma_p_vv_vari(a.vi_, b.vi_));
    }

    inline var gamma_p(const var& a,
                       double b) {
      return var(new gamma_p_vd_vari(a.vi_, b));
    }

    inline var gamma_p(double a,
                       const var& b) {
      return var(new gamma_p_dv_vari(a, b.vi_));
    }

  }
}
#endif
