#ifndef STAN_MATH_REV_FUN_LOG1P_EXP_HPP
#define STAN_MATH_REV_FUN_LOG1P_EXP_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/fun/log1p_exp.hpp>

namespace stan {
namespace math {

/**
 * Return the log of 1 plus the exponential of the specified
 * variable.
 * @tparam T Arithmetic or a type inheriting from `EigenBase`.
 * @param a The variable.
 */
template <typename T, require_stan_scalar_or_eigen_t<T>* = nullptr>
inline auto log1p_exp(const var_value<T>& a) {
  auto precomp_inv_logit = to_arena(as_array_or_scalar(inv_logit(a.val())));
  return make_callback_var(
      log1p_exp(a.val()), [a, precomp_inv_logit](auto& vi) mutable {
        as_array_or_scalar(a.adj())
            += as_array_or_scalar(vi.adj()) * precomp_inv_logit;
      });
}

}  // namespace math
}  // namespace stan
#endif
