#ifndef STAN_MATH_REV_FUN_ATAN_HPP
#define STAN_MATH_REV_FUN_ATAN_HPP

#include <stan/math/prim/fun/abs.hpp>
#include <stan/math/prim/fun/isinf.hpp>
#include <stan/math/prim/fun/isnan.hpp>
#include <stan/math/prim/fun/atan.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/hypot.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/abs.hpp>
#include <stan/math/rev/fun/cos.hpp>
#include <stan/math/rev/fun/atanh.hpp>
#include <stan/math/rev/fun/log.hpp>
#include <stan/math/rev/fun/sqrt.hpp>
#include <stan/math/rev/fun/is_inf.hpp>
#include <stan/math/rev/fun/is_nan.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

namespace internal {
class atan_vari : public op_v_vari {
 public:
  explicit atan_vari(vari* avi) : op_v_vari(std::atan(avi->val_), avi) {}
  void chain() { avi_->adj_ += adj_ / (1.0 + (avi_->val_ * avi_->val_)); }
};
}  // namespace internal

/**
 * Return the principal value of the arc tangent, in radians, of the
 * specified variable (cmath).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \arctan x = \frac{1}{1 + x^2}\f$.
 *
 *
   \f[
   \mbox{atan}(x) =
   \begin{cases}
     \arctan(x) & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{atan}(x)}{\partial x} =
   \begin{cases}
     \frac{\partial\, \arctan(x)}{\partial x} & \mbox{if } -\infty\leq x\leq
 \infty \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN} \end{cases} \f]

   \f[
   \frac{\partial \, \arctan(x)}{\partial x} = \frac{1}{x^2+1}
   \f]
 *
 * @param a Variable in range [-1, 1].
 * @return Arc tangent of variable, in radians.
 */
inline var atan(const var& a) { return var(new internal::atan_vari(a.vi_)); }

/**
 * Return the arc tangent of the complex argument.
 *
 * @param[in] z argument
 * @return arc tangent of the argument
 */
inline std::complex<var> atan(const std::complex<var>& z) {
  return stan::math::internal::complex_atan(z);
}

}  // namespace math
}  // namespace stan
#endif
