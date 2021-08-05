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
 * @param a Argument variable.
 * @return Inverse logit of argument.
 */
inline var inv_logit(const var& a) {
  return make_callback_var(inv_logit(a.val()), [a](auto& vi) mutable {
    a.adj() += vi.adj() * vi.val() * (1.0 - vi.val());
  });
}

/**
 * Return the inverse logit of the elements of x
 *
 * @tparam T type of x
 * @param x argument
 * @return elementwise inverse logit of x
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto inv_logit(const T& x) {
  return make_callback_var(inv_logit(x.val()), [x](const auto& vi) mutable {
    x.adj() += (vi.adj().array() * vi.val().array() * (1.0 - vi.val().array()))
                   .matrix();
  });
}

}  // namespace math
}  // namespace stan
#endif
