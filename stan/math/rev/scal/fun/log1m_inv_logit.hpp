#ifndef STAN_MATH_REV_SCAL_FUN_LOG1M_INV_LOGIT_HPP
#define STAN_MATH_REV_SCAL_FUN_LOG1M_INV_LOGIT_HPP

#include <stan/math/rev/core/precomp_v_vari.hpp>
#include <stan/math/rev/scal/fun/log1p.hpp>
#include <stan/math/rev/scal/fun/exp.hpp>

namespace stan {
  namespace math {

    /**
     * Return the natural logarithm of one minus the inverse logit of
     * the specified argument.
     *
     * @param u argument
     * @return log of one minus the inverse logit of the argument
     */
    inline var log1m_inv_logit(const var& u) {
      double val = log1m_inv_logit(u.val());
      return var(new precomp_v_vari(val, u.vi_, -inv_logit(u.val())));
    }

  }
}
#endif
