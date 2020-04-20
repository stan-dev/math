#ifndef STAN_MATH_REV_FUN_TAN_HPP
#define STAN_MATH_REV_FUN_TAN_HPP

#include <stan/math/prim/fun/isinf.hpp>
#include <stan/math/prim/fun/tan.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/is_inf.hpp>
#include <stan/math/rev/fun/sinh.hpp>
#include <stan/math/rev/fun/tanh.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

namespace internal {
class tan_vari : public op_v_vari {
 public:
  explicit tan_vari(vari* avi) : op_v_vari(std::tan(avi->val_), avi) {}
  void chain() { avi_->adj_ += adj_ * (1.0 + val_ * val_); }
};
}  // namespace internal

/**
 * Return the tangent of a radian-scaled variable (cmath).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \tan x = \sec^2 x\f$.
 *
 *
   \f[
   \mbox{tan}(x) =
   \begin{cases}
     \tan(x) & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{tan}(x)}{\partial x} =
   \begin{cases}
     \sec^2(x) & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Variable for radians of angle.
 * @return Tangent of variable.
 */
inline var tan(const var& a) { return var(new internal::tan_vari(a.vi_)); }

/**
 * Return the tangent of the complex argument.
 *
 * @param[in] z argument
 * @return tangent of the argument
 */
inline std::complex<var> tan(const std::complex<var>& z) {
  return stan::math::internal::complex_tan(z);
}

}  // namespace math
}  // namespace stan
#endif
