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
template <typename T, require_stan_scalar_or_eigen_t<T>* = nullptr>
inline auto inv_logit(const var_value<T>& a) {
  return make_callback_var(inv_logit(a.val()), [a](auto& vi) mutable {
    as_array_or_scalar(a).adj() += as_array_or_scalar(vi.adj())
                                   * as_array_or_scalar(vi.val())
                                   * (1.0 - as_array_or_scalar(vi.val()));
  });
}

}  // namespace math
}  // namespace stan
#endif
