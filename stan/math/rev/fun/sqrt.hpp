#ifndef STAN_MATH_REV_FUN_SQRT_HPP
#define STAN_MATH_REV_FUN_SQRT_HPP

#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/atan2.hpp>
#include <stan/math/rev/fun/cos.hpp>
#include <stan/math/rev/fun/hypot.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

namespace internal {
class sqrt_vari : public op_v_vari {
 public:
  explicit sqrt_vari(vari* avi) : op_v_vari(std::sqrt(avi->val_), avi) {}
  void chain() { avi_->adj_ += adj_ / (2.0 * val_); }
};
}  // namespace internal

/**
 * Return the square root of the specified variable (cmath).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \sqrt{x} = \frac{1}{2 \sqrt{x}}\f$.
 *
   \f[
   \mbox{sqrt}(x) =
   \begin{cases}
     \textrm{NaN} & x < 0 \\
     \sqrt{x} & \mbox{if } x\geq 0\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{sqrt}(x)}{\partial x} =
   \begin{cases}
     \textrm{NaN} & x < 0 \\
     \frac{1}{2\sqrt{x}} & x\geq 0\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Variable whose square root is taken.
 * @return Square root of variable.
 */
inline var sqrt(const var& a) { return var(new internal::sqrt_vari(a.vi_)); }

/**
 * Return the square root of the complex argument.
 *
 * @param[in] z argument
 * @return square root of the argument
 */
inline std::complex<var> sqrt(const std::complex<var>& z) {
  return internal::complex_sqrt(z);
}

}  // namespace math
}  // namespace stan
#endif
