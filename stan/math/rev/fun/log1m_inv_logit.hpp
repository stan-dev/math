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
 * @tparam T Arithmetic or a type inheriting from `EigenBase`.
 * @param u argument
 * @return log of one minus the inverse logit of the argument
 */
inline auto log1m_inv_logit(const var u) {
  auto precomp_inv_logit = -inv_logit(u.val());
  return make_callback_var(log1m_inv_logit(u.val()),
                           [u, precomp_inv_logit](auto& vi) mutable {
                             u.adj() += vi.adj() * precomp_inv_logit;
                           });
}

template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto log1m_inv_logit(const T& u) {
  auto u_arena = to_arena(u);
  auto precomp_inv_logit = to_arena(-inv_logit(u_arena.val()).array());
  return make_callback_rev_matrix<T>(
      log1m_inv_logit(u_arena.val()),
      [u_arena, precomp_inv_logit](auto&& vi) mutable {
        u_arena.adj().array() += vi.adj().array() * precomp_inv_logit;
      });
}

}  // namespace math
}  // namespace stan
#endif
