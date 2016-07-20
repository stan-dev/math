#ifndef STAN_MATH_REV_SCAL_FUN_SQUARED_DISTANCE_HPP
#define STAN_MATH_REV_SCAL_FUN_SQUARED_DISTANCE_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/scal/fun/squared_distance.hpp>

namespace stan {
  namespace math {

    class scal_squared_distance_vv_vari : public op_vv_vari {
    public:
      scal_squared_distance_vv_vari(vari* avi, vari* bvi) :
        op_vv_vari(stan::math::squared_distance(avi->val_, bvi->val_),
                   avi, bvi) {
      }
      void chain() {
        double diff = avi_->val_ - bvi_->val_;
        avi_->adj_ += adj_ * 2.0 * diff;
        bvi_->adj_ -= adj_ * 2.0 * diff;
      }
    };
    class scal_squared_distance_vd_vari : public op_vd_vari {
    public:
      scal_squared_distance_vd_vari(vari* avi, double b) :
        op_vd_vari(stan::math::squared_distance(avi->val_, b),
                   avi, b) {
      }
      void chain() {
        avi_->adj_ += adj_ * 2 * (avi_->val_ - bd_);
      }
    };
    class scal_squared_distance_dv_vari : public op_dv_vari {
    public:
      scal_squared_distance_dv_vari(double a, vari* bvi) :
        op_dv_vari(stan::math::squared_distance(a, bvi->val_),
                   a, bvi) {
      }
      void chain() {
        bvi_->adj_ -= adj_ * 2 * (ad_ - bvi_->val_);
      }
    };

    /**
     * Returns the log sum of exponentials.
     */
    inline var squared_distance(const stan::math::var& a,
                            const stan::math::var& b) {
      return var(new scal_squared_distance_vv_vari(a.vi_, b.vi_));
    }
    /**
     * Returns the log sum of exponentials.
     */
    inline var squared_distance(const stan::math::var& a,
                            const double& b) {
      return var(new scal_squared_distance_vd_vari(a.vi_, b));
    }
    /**
     * Returns the log sum of exponentials.
     */
    inline var squared_distance(const double& a,
                                const stan::math::var& b) {
      return var(new scal_squared_distance_dv_vari(a, b.vi_));
    }
  }
}
#endif
