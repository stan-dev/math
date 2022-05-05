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
inline auto log1p_exp(const var a) {
  auto precomp_inv_logit = inv_logit(a.val());
  return make_callback_var(log1p_exp(a.val()),
                           [a, precomp_inv_logit](auto& vi) mutable {
                             a.adj() += vi.adj() * precomp_inv_logit;
                           });
}

template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto log1p_exp(const T& a) {
  auto a_arena = to_arena(a);
  auto precomp_inv_logit = to_arena(inv_logit(a_arena.val()).array());
  return make_callback_rev_matrix<T>(
      log1p_exp(a_arena.val()),
      [a_arena, precomp_inv_logit](auto&& vi) mutable {
        a_arena.adj().array() += vi.adj().array() * precomp_inv_logit;
      });
}

}  // namespace math
}  // namespace stan
#endif
