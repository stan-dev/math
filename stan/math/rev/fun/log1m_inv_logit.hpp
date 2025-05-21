#ifndef STAN_MATH_REV_FUN_LOG1M_INV_LOGIT_HPP
#define STAN_MATH_REV_FUN_LOG1M_INV_LOGIT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/inv_logit.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/log1m_inv_logit.hpp>

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
template <typename T, require_stan_scalar_or_eigen_t<T>* = nullptr>
inline auto log1m_inv_logit(const var_value<T>& u) {
  auto precomp_inv_logit = to_arena(as_array_or_scalar(-inv_logit(u.val())));
  return make_callback_var(
      log1m_inv_logit(u.val()), [u, precomp_inv_logit](auto& vi) mutable {
        as_array_or_scalar(u.adj())
            += as_array_or_scalar(vi.adj()) * precomp_inv_logit;
      });
}

}  // namespace math
}  // namespace stan
#endif
