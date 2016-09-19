#ifndef STAN_MATH_REV_SCAL_FUN_FREXP_HPP
#define STAN_MATH_REV_SCAL_FUN_FREXP_HPP

#include <stan/math/rev/core.hpp>
#include <cmath>
#include <limits>

namespace stan {
  namespace math {

	/**
	 * Variable implementation of frexp. The derivative of frexp
	 * is 0, because frexp returns a discontinuous object. 
	 */
      class frexp_vari : public op_v_vari {
      public:
      /**
       * Construct the variable implementation of frexp.
       */
        explicit frexp_vari(vari* avi, int* b) :
        op_v_vari(std::frexp(avi->val_, b), avi) {
        }
      /**
       * Assign a quiet_NaN value to the adjoint, because frexp
       * returns a discontinuous object.
       */
        void chain() {
          if (unlikely(boost::math::isnan(avi_->val_)))
            avi_->adj_ = std::numeric_limits<double>::quiet_NaN();
        }
      };


    /**
     * Decomposes a float-like variable into a normalized 
     * fraction and an integral power of two. (cmath)
     *
     * @param[in] a The variable that will be decomposed.
     * @param[out] b Pointer to an integer that will stored the
     * integral power of two.
     * @return Normalized fraction
     */
      inline var frexp(const var& a, int* b) {
        return var(new frexp_vari(a.vi_, b));
      }

  }
}

#endif
