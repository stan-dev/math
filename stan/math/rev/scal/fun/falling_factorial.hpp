#ifndef STAN_MATH_REV_SCAL_FUN_FALLING_FACTORIAL_HPP
#define STAN_MATH_REV_SCAL_FUN_FALLING_FACTORIAL_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>
#include <stan/math/prim/scal/fun/falling_factorial.hpp>

namespace stan {
  namespace math {

    namespace {

      class falling_factorial_vv_vari : public op_vv_vari {
      public:
        falling_factorial_vv_vari(vari* avi, vari* bvi) :
          op_vv_vari(falling_factorial(avi->val_, bvi->val_),
                     avi, bvi) {
        }
        void chain() {
          avi_->adj_ += adj_
            * val_
            * (digamma(avi_->val_ + 1)
               - digamma(avi_->val_ - bvi_->val_ + 1));
          bvi_->adj_ += adj_
            * val_
            * digamma(avi_->val_ - bvi_->val_ + 1);
        }
      };

      class falling_factorial_vd_vari : public op_vd_vari {
      public:
        falling_factorial_vd_vari(vari* avi, double b) :
          op_vd_vari(falling_factorial(avi->val_, b), avi, b) {
        }
        void chain() {
          avi_->adj_ += adj_
            * val_
            * (digamma(avi_->val_ + 1)
               - digamma(avi_->val_ - bd_ + 1));
        }
      };

      class falling_factorial_dv_vari : public op_dv_vari {
      public:
        falling_factorial_dv_vari(double a, vari* bvi) :
          op_dv_vari(falling_factorial(a, bvi->val_), a, bvi) {
        }
        void chain() {
          bvi_->adj_ += adj_
            * val_
            * digamma(ad_ - bvi_->val_ + 1);
        }
      };
    }

    inline var falling_factorial(const var& a,
                                 double b) {
      return var(new falling_factorial_vd_vari(a.vi_, b));
    }

    inline var falling_factorial(const var& a,
                                 const var& b) {
      return var(new falling_factorial_vv_vari(a.vi_, b.vi_));
    }

    inline var falling_factorial(double a,
                                 const var& b) {
      return var(new falling_factorial_dv_vari(a, b.vi_));
    }

  }
}
#endif
