#ifndef STAN_MATH_REV_FUN_ACOSH_HPP
#define STAN_MATH_REV_FUN_ACOSH_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/abs.hpp>
#include <stan/math/rev/fun/arg.hpp>
#include <stan/math/rev/fun/cosh.hpp>
#include <stan/math/rev/fun/is_nan.hpp>
#include <stan/math/rev/fun/log.hpp>
#include <stan/math/rev/fun/polar.hpp>
#include <stan/math/rev/fun/sqrt.hpp>
#include <stan/math/rev/fun/value_of_rec.hpp>
#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/isnan.hpp>
#include <stan/math/prim/fun/acosh.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * The inverse hyperbolic cosine function for variables (C99).
 *
 * For non-variable function, see ::acosh().
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \mbox{acosh}(x) = \frac{x}{x^2 - 1}\f$.
 *
 *
   \f[
   \mbox{acosh}(x) =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < 1 \\
     \cosh^{-1}(x) & \mbox{if } x \geq 1 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{acosh}(x)}{\partial x} =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < 1 \\
     \frac{\partial\, \cosh^{-1}(x)}{\partial x} & \mbox{if } x \geq 1 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \cosh^{-1}(x)=\ln\left(x+\sqrt{x^2-1}\right)
   \f]

   \f[
   \frac{\partial \, \cosh^{-1}(x)}{\partial x} = \frac{1}{\sqrt{x^2-1}}
   \f]
 *
 * @param x The variable.
 * @return Inverse hyperbolic cosine of the variable.
 */
inline var acosh(const var& x) {
  return make_callback_var(acosh(x.val()), [x](const auto& vi) mutable {
    x.adj() += vi.adj() / std::sqrt(x.val() * x.val() - 1.0);
  });
}
/**
 * The inverse hyperbolic cosine function for variables (C99).
 *
 * For non-variable function, see ::acosh().
 *
 * @tparam Varmat a `var_value` with inner Eigen type
 * @param x The variable
 * @return Inverse hyperbolic cosine of the variable.
 */
template <typename VarMat, require_var_matrix_t<VarMat>* = nullptr>
inline auto acosh(const VarMat& x) {
  return make_callback_var(
      x.val().unaryExpr([](const auto x) { return acosh(x); }),
      [x](const auto& vi) mutable {
        x.adj().array()
            += vi.adj().array() / (x.val().array().square() - 1.0).sqrt();
      });
}

/**
 * Return the hyperbolic arc cosine of the complex argument.
 *
 * @param[in] z argument
 * @return hyperbolic arc cosine of the argument
 */
inline std::complex<var> acosh(const std::complex<var>& z) {
  return stan::math::internal::complex_acosh(z);
}

}  // namespace math
}  // namespace stan
#endif
