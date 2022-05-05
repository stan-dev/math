#ifndef STAN_MATH_REV_FUN_INV_LOGIT_HPP
#define STAN_MATH_REV_FUN_INV_LOGIT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>

namespace stan {
namespace math {

/**
 * The inverse logit function for variables (stan).
 *
 * See inv_logit() for the double-based version.
 *
 * The derivative of inverse logit is
 *
 * \f$\frac{d}{dx} \mbox{logit}^{-1}(x) = \mbox{logit}^{-1}(x) (1 -
 * \mbox{logit}^{-1}(x))\f$.
 *
 * @tparam T Arithmetic or a type inheriting from `EigenBase`.
 * @param a Argument variable.
 * @return Inverse logit of argument.
 */
inline auto inv_logit(const var a) {
  return make_callback_var(inv_logit(a.val()), [a](auto& vi) mutable {
    a.adj() += vi.adj() * vi.val() * (1.0 - vi.val());
  });
}

template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto inv_logit(const T& x) {
  auto x_arena = to_arena(x);
  return make_callback_rev_matrix<T>(inv_logit(x_arena.val()), [x_arena](auto& vi) mutable {
    x_arena.adj().array() += vi.adj().array() * vi.val().array() * (1.0 - vi.val().array());
  });
}

}  // namespace math
}  // namespace stan
#endif
