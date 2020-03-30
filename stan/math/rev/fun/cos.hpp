#ifndef STAN_MATH_REV_FUN_COS_HPP
#define STAN_MATH_REV_FUN_COS_HPP

#include <stan/math/prim/fun/abs.hpp>
#include <stan/math/prim/fun/cos.hpp>
#include <stan/math/prim/fun/isinf.hpp>
#include <stan/math/prim/fun/isfinite.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/abs.hpp>
#include <stan/math/rev/fun/cosh.hpp>
#include <stan/math/rev/fun/sinh.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

namespace internal {
class cos_vari : public op_v_vari {
 public:
  explicit cos_vari(vari* avi) : op_v_vari(std::cos(avi->val_), avi) {}
  void chain() { avi_->adj_ -= adj_ * std::sin(avi_->val_); }
};
}  // namespace internal

/**
 * Return the cosine of a radian-scaled variable (cmath).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \cos x = - \sin x\f$.
 *
 *
   \f[
   \mbox{cos}(x) =
   \begin{cases}
     \cos(x) & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{cos}(x)}{\partial x} =
   \begin{cases}
     -\sin(x) & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Variable for radians of angle.
 * @return Cosine of variable.
 */
inline var cos(const var& a) { return var(new internal::cos_vari(a.vi_)); }

/**
 * Return the cosine of the complex argument.
 *
 * @param[in] z argument
 * @return cosine of the argument
 */
inline std::complex<var> cos(const std::complex<var>& z) {
  return stan::math::internal::complex_cos(z);
}

}  // namespace math
}  // namespace stan
#endif
