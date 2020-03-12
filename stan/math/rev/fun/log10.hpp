#ifndef STAN_MATH_REV_FUN_LOG10_HPP
#define STAN_MATH_REV_FUN_LOG10_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log10.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/atan2.hpp>
#include <stan/math/rev/fun/hypot.hpp>
#include <stan/math/rev/fun/log.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

namespace internal {
class log10_vari : public op_v_vari {
 public:
  const double exp_val_;
  explicit log10_vari(vari* avi)
      : op_v_vari(std::log10(avi->val_), avi), exp_val_(avi->val_) {}
  void chain() { avi_->adj_ += adj_ / (LOG_TEN * exp_val_); }
};
}  // namespace internal

/**
 * Return the base 10 log of the specified variable (cmath).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \log_{10} x = \frac{1}{x \log 10}\f$.
 *
 *
   \f[
   \mbox{log10}(x) =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < 0\\
     \log_{10}(x) & \mbox{if } x \geq 0 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log10}(x)}{\partial x} =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < 0\\
     \frac{1}{x \ln10} & \mbox{if } x\geq 0 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Variable whose log is taken.
 * @return Base 10 log of variable.
 */
inline var log10(const var& a) { return var(new internal::log10_vari(a.vi_)); }

/**
 * Return the base 10 logarithm of the specified complex number.
 *
 * @param z complex argument
 * @return base 10 log of argument
 */
inline std::complex<var> log10(const std::complex<var>& z) {
  return internal::complex_log10(z);
}

}  // namespace math
}  // namespace stan
#endif
