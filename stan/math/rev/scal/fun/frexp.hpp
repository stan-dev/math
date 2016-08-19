#ifndef STAN_MATH_REV_SCAL_FUN_FREXP_HPP
#define STAN_MATH_REV_SCAL_FUN_FREXP_HPP

#include <stan/math/rev/core.hpp>
#include <cmath>

namespace stan {
  namespace math {
      
      class frexp_vari : public op_v_vari {
      public:
          explicit frexp_vari(vari* avi, int* b) :
          op_v_vari(::frexp(avi->val_, b), avi) {
          }
          void chain() {
              if (unlikely(boost::math::isnan(avi_->val_)))
                  avi_->adj_ = std::numeric_limits<double>::quiet_NaN();
          }
      };


    /**
     * Decomposes a float-like variable into a normalized 
     * fraction and an integral power of two.
     *
     * @param the variable that will be decomposed.
     * @param pointer to an integer that will stored the
     *        integral power of two.
     * @return normalized fraction
     */
      inline var frexp(const var& a, int* b) {
          return var(new frexp_vari(a.vi_, b));
      }

  }
}

#endif
