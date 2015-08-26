#ifndef STAN_MATH_REV_SCAL_FUN_INV_PHI_HPP
#define STAN_MATH_REV_SCAL_FUN_INV_PHI_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/scal/fun/inv_Phi.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>

namespace stan {
  namespace math {

    namespace {
      class inv_Phi_vari : public op_v_vari {
      private:
        double m_z;
      public:
        explicit inv_Phi_vari(vari* avi) :
          op_v_vari(m_z = stan::math::inv_Phi(avi->val_), avi) {
        }
        void chain() {
          static const double NEG_HALF = -0.5;
          avi_->adj_ += adj_
            * SQRT_2_TIMES_SQRT_PI
            / std::exp(NEG_HALF * m_z * m_z);
        }
      };
    }

    /**
     * The inverse of unit normal cumulative density function.
     *
     * See stan::math::inv_Phi() for the double-based version.
     *
     * The derivative is the reciprocal of unit normal density function,
     *
     * @param p Probability
     * @return The unit normal inverse cdf evaluated at p
     */
    inline var inv_Phi(const stan::math::var& p) {
      return var(new inv_Phi_vari(p.vi_));
    }

  }
}
#endif
