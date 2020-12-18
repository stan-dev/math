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
 * @param x Variable in range [-1, 1].
 * @return Arc tangent of variable, in radians.
 */
inline var atan(const var& x) {
  return make_callback_var(std::atan(x.val()), [x](const auto& vi) mutable {
    x.adj() += vi.adj() / (1.0 + (x.val() * x.val()));
  });
}

/**
 * Return the principal value of the arc tangent, in radians, of the
 * specified variable (cmath).
 *
 *
 * @tparam Varmat a `var_value` with inner Eigen type
 * @param x Variable in range [-1, 1].
 * @return Arc tangent of variable, in radians.
 */
template <typename VarMat, require_var_matrix_t<VarMat>* = nullptr>
inline auto atan(const VarMat& x) {
  return make_callback_var(
      x.val().array().atan().matrix(), [x](const auto& vi) mutable {
        x.adj().array()
            += vi.adj().array() / (1.0 + (x.val().array().square()));
      });
}
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
