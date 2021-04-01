#ifndef STAN_MATH_REV_FUN_LOG_INV_LOGIT_HPP
#define STAN_MATH_REV_FUN_LOG_INV_LOGIT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/fun/log_inv_logit.hpp>

namespace stan {
namespace math {

/**
 * Return the natural logarithm of the inverse logit of the
 * specified argument.
 *
 * @param u argument
 * @return log inverse logit of the argument
 */
inline var log_inv_logit(const var& u) {
  return make_callback_var(log_inv_logit(u.val()), [u](auto& vi) mutable {
    u.adj() += vi.adj() * inv_logit(-u.val());
  });
}

}  // namespace math
}  // namespace stan
#endif
