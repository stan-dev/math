#ifndef STAN_MATH_REV_FUN_PHI_HPP
#define STAN_MATH_REV_FUN_PHI_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/Phi.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * The unit normal cumulative density function for variables (stan).
 *
 * See Phi() for the double-based version.
 *
 * The derivative is the unit normal density function,
 *
 * \f$\frac{d}{dx} \Phi(x) = \mbox{\sf Norm}(x|0, 1) = \frac{1}{\sqrt{2\pi}}
 \exp(-\frac{1}{2} x^2)\f$.
 *
 *
   \f[
   \mbox{Phi}(x) =
   \begin{cases}
     0 & \mbox{if } x < -37.5 \\
     \Phi(x) & \mbox{if } -37.5 \leq x \leq 8.25 \\
     1 & \mbox{if } x > 8.25 \\[6pt]
     \textrm{error} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{Phi}(x)}{\partial x} =
   \begin{cases}
     0 & \mbox{if } x < -27.5 \\
     \frac{\partial\, \Phi(x)}{\partial x} & \mbox{if } -27.5 \leq x \leq 27.5
 \\ 0 & \mbox{if } x > 27.5 \\[6pt] \textrm{error} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \Phi(x) = \frac{1}{\sqrt{2\pi}} \int_{0}^{x} e^{-t^2/2} dt
   \f]

   \f[
   \frac{\partial \, \Phi(x)}{\partial x} = \frac{e^{-x^2/2}}{\sqrt{2\pi}}
   \f]
 *
 * @param a Variable argument.
 * @return The unit normal cdf evaluated at the specified argument.
 */
inline var Phi(const var& a) {
  return make_callback_var(Phi(a.val()), [a](auto& vi) mutable {
    a.adj() += vi.adj() * INV_SQRT_TWO_PI * std::exp(-0.5 * a.val() * a.val());
  });
}

/**
 * Elementwise unit normal cumulative density function for varmat types.
 *
 * @tparam T a `var_value` with inner Eigen type
 * @param a input
 * @return The vectorized unit normal cdf
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto Phi(const T& a) {
  return make_callback_var(Phi(a.val()), [a](auto& vi) mutable {
    a.adj().array() += vi.adj().array() * INV_SQRT_TWO_PI
                       * (-0.5 * a.val().array().square()).exp();
  });
}

}  // namespace math
}  // namespace stan
#endif
