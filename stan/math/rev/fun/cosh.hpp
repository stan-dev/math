#ifndef STAN_MATH_REV_FUN_COSH_HPP
#define STAN_MATH_REV_FUN_COSH_HPP

#include <stan/math/prim/fun/cosh.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/operator_multiplication.hpp>
#include <stan/math/rev/fun/exp.hpp>
#include <stan/math/rev/fun/sin.hpp>
#include <stan/math/rev/fun/sinh.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

namespace internal {
class cosh_vari : public op_v_vari {
 public:
  explicit cosh_vari(vari* avi) : op_v_vari(std::cosh(avi->val_), avi) {}
  void chain() { avi_->adj_ += adj_ * std::sinh(avi_->val_); }
};
}  // namespace internal

/**
 * Return the hyperbolic cosine of the specified variable (cmath).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \cosh x = \sinh x\f$.
 *
 *
   \f[
   \mbox{cosh}(x) =
   \begin{cases}
     \cosh(x) & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{cosh}(x)}{\partial x} =
   \begin{cases}
     \sinh(x) & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Variable.
 * @return Hyperbolic cosine of variable.
 */
inline var cosh(const var& a) { return var(new internal::cosh_vari(a.vi_)); }

/**
 * Return the hyperbolic cosine of the complex argument.
 *
 * @param[in] z argument
 * @return hyperbolic cosine of the argument
 */
inline std::complex<var> cosh(const std::complex<var>& z) {
  return stan::math::internal::complex_cosh(z);
}

}  // namespace math
}  // namespace stan
#endif
