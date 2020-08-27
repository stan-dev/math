#ifndef STAN_MATH_REV_FUN_LOG_HPP
#define STAN_MATH_REV_FUN_LOG_HPP

#include <stan/math/prim/fun/log.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/abs.hpp>
#include <stan/math/rev/fun/arg.hpp>
#include <stan/math/rev/fun/atan2.hpp>
#include <stan/math/rev/fun/cos.hpp>
#include <stan/math/rev/fun/is_inf.hpp>
#include <stan/math/rev/fun/is_nan.hpp>
#include <stan/math/rev/fun/sqrt.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {
class log_vari final : public op_v_vari {
 public:
  explicit log_vari(vari* avi) : op_v_vari(std::log(avi->val_), avi) {}
  void chain() { avi_->adj_ += adj_ / avi_->val_; }
};
}  // namespace internal

/**
 * Return the natural log of the specified variable (cmath).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \log x = \frac{1}{x}\f$.
 *
   \f[
   \mbox{log}(x) =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < 0\\
     \ln(x) & \mbox{if } x \geq 0 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log}(x)}{\partial x} =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < 0\\
     \frac{1}{x} & \mbox{if } x\geq 0 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Variable whose log is taken.
 * @return Natural log of variable.
 */
inline var log(const var& a) { return var(new internal::log_vari(a.vi_)); }

/**
 * Return the natural logarithm (base e) of the specified complex argument.
 *
 * @param z complex argument
 * @return natural logarithm of argument
 */
inline std::complex<var> log(const std::complex<var>& z) {
  return internal::complex_log(z);
}

}  // namespace math
}  // namespace stan
#endif
