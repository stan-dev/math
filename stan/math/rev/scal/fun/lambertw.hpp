#ifndef STAN_MATH_REV_SCAL_FUN_LAMBERT_W_HPP
#define STAN_MATH_REV_SCAL_FUN_LAMBERT_W_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/scal/fun/lambertw.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

namespace internal {
class lambert_w0_vari : public op_v_vari {
 public:
  explicit lambert_w0_vari(vari* avi) : op_v_vari(lambert_w0(avi->val_), avi) {}
  void chain() { avi_->adj_ += adj_ * (this->val_ / (avi->val_ * (1 + this->val_))); }
};

class lambert_w1_vari : public op_v_vari {
 public:
  explicit lambert_w1_vari(vari* avi) : op_v_vari(lambert_w1(avi->val_), avi) {}
  void chain() { avi_->adj_ += adj_ * (this->val_ / (avi->val_ * (1 + this->val_))); }
};

}  // namespace internal

/**
 * The exponentiation of the specified variable minus 1 (C99).
 *
 * The derivative is given by
 *
 * \f$\frac{d}{dx} \exp(a) - 1 = \exp(a)\f$.
 *
 *
   \f[
   \mbox{lambert_w0}(x) =
   \begin{cases}
     e^x-1 & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{lambert_w0}(x)}{\partial x} =
   \begin{cases}
     e^x & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a The variable.
 * @return Two to the power of the specified variable.
 */
inline var lambert_w0(const var& a) { return var(new internal::lambert_w0_vari(a.vi_)); }
inline var lambert_w1(const var& a) { return var(new internal::lambert_w1_vari(a.vi_)); }

}  // namespace math
}  // namespace stan
#endif
