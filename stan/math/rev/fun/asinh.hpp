#ifndef STAN_MATH_REV_FUN_ASINH_HPP
#define STAN_MATH_REV_FUN_ASINH_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/asinh.hpp>
#include <stan/math/prim/fun/isinf.hpp>
#include <stan/math/prim/fun/is_inf.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/value_of_rec.hpp>
#include <stan/math/rev/fun/abs.hpp>
#include <stan/math/rev/fun/arg.hpp>
#include <stan/math/rev/fun/cosh.hpp>
#include <stan/math/rev/fun/is_inf.hpp>
#include <stan/math/rev/fun/is_nan.hpp>
#include <stan/math/rev/fun/log.hpp>
#include <stan/math/rev/fun/polar.hpp>
#include <stan/math/rev/fun/sqrt.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

namespace internal {
class asinh_vari : public op_v_vari {
 public:
  asinh_vari(double val, vari* avi) : op_v_vari(val, avi) {}
  void chain() {
    avi_->adj_ += adj_ / std::sqrt(avi_->val_ * avi_->val_ + 1.0);
  }
};
}  // namespace internal

/**
 * The inverse hyperbolic sine function for variables (C99).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \mbox{asinh}(x) = \frac{x}{x^2 + 1}\f$.
 *
 *
   \f[
   \mbox{asinh}(x) =
   \begin{cases}
     \sinh^{-1}(x) & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{asinh}(x)}{\partial x} =
   \begin{cases}
     \frac{\partial\, \sinh^{-1}(x)}{\partial x} & \mbox{if } -\infty\leq x\leq
 \infty \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN} \end{cases} \f]

   \f[
   \sinh^{-1}(x)=\ln\left(x+\sqrt{x^2+1}\right)
   \f]

   \f[
   \frac{\partial \, \sinh^{-1}(x)}{\partial x} = \frac{1}{\sqrt{x^2+1}}
   \f]
 *
 * @param a The variable.
 * @return Inverse hyperbolic sine of the variable.
 */
inline var asinh(const var& a) {
  return var(new internal::asinh_vari(asinh(a.val()), a.vi_));
}

/**
 * Return the hyperbolic arcsine of the complex argument.
 *
 * @param[in] z argument
 * @return hyperbolic arcsine of the argument
 */
inline std::complex<var> asinh(const std::complex<var>& z) {
  return stan::math::internal::complex_asinh(z);
}

}  // namespace math
}  // namespace stan
#endif
