#ifndef STAN_MATH_REV_FUN_ASINH_HPP
#define STAN_MATH_REV_FUN_ASINH_HPP

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
#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/asinh.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

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
 * @param x The variable.
 * @return Inverse hyperbolic sine of the variable.
 */
inline var asinh(const var& x) {
  return make_callback_var(std::asinh(x.val()), [x](const auto& vi) mutable {
    x.adj() += vi.adj() / std::sqrt(x.val() * x.val() + 1.0);
  });
}

/**
 * The inverse hyperbolic sine function for variables (C99).
 *
 * @tparam Varmat a `var_value` with inner Eigen type
 * @param x The variable.
 * @return Inverse hyperbolic sine of the variable.
 */
template <typename VarMat, require_var_matrix_t<VarMat>* = nullptr>
inline auto asinh(const VarMat& x) {
  return make_callback_var(
      x.val().unaryExpr([](const auto x) { return asinh(x); }),
      [x](const auto& vi) mutable {
        x.adj().array()
            += vi.adj().array() / (x.val().array().square() + 1.0).sqrt();
      });
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
