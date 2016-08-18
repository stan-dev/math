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
     * Return the sine of a radian-scaled variable (cmath).
     *
     * The derivative is defined by
     *
     * \f$\frac{d}{dx} \sin x = \cos x\f$.
     *
     *
       \f[
       \mbox{sin}(x) =
       \begin{cases}
         \sin(x) & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
         \textrm{NaN} & \mbox{if } x = \textrm{NaN}
       \end{cases}
       \f]

       \f[
       \frac{\partial\, \mbox{sin}(x)}{\partial x} =
       \begin{cases}
         \cos(x) & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
         \textrm{NaN} & \mbox{if } x = \textrm{NaN}
       \end{cases}
       \f]
     *
     * @param a Variable for radians of angle.
     * @return Sine of variable.
     */
      inline var frexp(const var& a, int* b) {
          return var(new frexp_vari(a.vi_, b));
      }

  }
}

#endif
