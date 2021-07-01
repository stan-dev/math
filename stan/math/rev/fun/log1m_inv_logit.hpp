#ifndef STAN_MATH_REV_FUN_LOG1M_INV_LOGIT_HPP
#define STAN_MATH_REV_FUN_LOG1M_INV_LOGIT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/fun/log1m_inv_logit.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/rev/core.hpp>

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
  auto precomp_inv_logit = -inv_logit(u.val());
  return make_callback_var(log1m_inv_logit(u.val()),
                           [u, precomp_inv_logit](auto& vi) mutable {
                             u.adj() += vi.adj() * precomp_inv_logit;
                           });
}

template <typename T, require_eigen_t<T>* = nullptr>
inline auto log1m_inv_logit(const var_value<T>& u) {
  auto precomp_inv_logit = to_arena(-inv_logit(u.val()).array());
  return make_callback_var(log1m_inv_logit(u.val()),
                           [u, precomp_inv_logit](auto& vi) mutable {
                             u.adj().array() += vi.adj().array() * precomp_inv_logit;
                           });
}

}  // namespace math
}  // namespace stan
#endif
