#ifndef STAN_MATH_REV_FUN_LOGIT_HPP
#define STAN_MATH_REV_FUN_LOGIT_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/fun/logit.hpp>

namespace stan {
namespace math {

/**
 * Return the log odds of the specified argument.
 *
 * @param u argument
 * @return log odds of argument
 */
template <typename T>
inline var_value<T> logit(const var_value<T>& u) {
  return {new precomp_v_vari<T, vari_value<T>>(
      logit(u.val()), u.vi_, 1.0 / (u.val() - u.val() * u.val()))};
}

}  // namespace math
}  // namespace stan
#endif
