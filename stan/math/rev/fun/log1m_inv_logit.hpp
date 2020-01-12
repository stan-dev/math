#ifndef STAN_MATH_REV_FUN_LOG1M_INV_LOGIT_HPP
#define STAN_MATH_REV_FUN_LOG1M_INV_LOGIT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/fun/log1m_inv_logit.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/precomp_v_vari.hpp>

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
  return var(
      new precomp_v_vari(log1m_inv_logit(u.val()), u.vi_, -inv_logit(u.val())));
}

}  // namespace math
}  // namespace stan
#endif
