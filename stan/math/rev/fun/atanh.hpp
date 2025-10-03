#ifndef STAN_MATH_REV_FUN_ATANH_HPP
#define STAN_MATH_REV_FUN_ATANH_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of_rec.hpp>
#include <stan/math/rev/fun/atan2.hpp>
#include <stan/math/rev/fun/cosh.hpp>
#include <stan/math/rev/fun/log.hpp>
#include <stan/math/rev/fun/sinh.hpp>
#include <stan/math/rev/fun/hypot.hpp>
#include <stan/math/prim/fun/atanh.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * The inverse hyperbolic tangent function for variables (C99).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \mbox{atanh}(x) = \frac{1}{1 - x^2}\f$.
 *
   \f[
   \mbox{atanh}(x) =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < -1\\
     \tanh^{-1}(x) & \mbox{if } -1\leq x \leq 1 \\
     \textrm{NaN} & \mbox{if } x > 1\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{atanh}(x)}{\partial x} =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < -1\\
     \frac{\partial\, \tanh^{-1}(x)}{\partial x} & \mbox{if } -1\leq x\leq 1 \\
     \textrm{NaN} & \mbox{if } x > 1\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \tanh^{-1}(x)=\frac{1}{2}\ln\left(\frac{1+x}{1-x}\right)
   \f]

   \f[
   \frac{\partial \, \tanh^{-1}(x)}{\partial x} = \frac{1}{1-x^2}
   \f]
   *
   * @param x The variable.
   * @return Inverse hyperbolic tangent of the variable.
   * @throw std::domain_error if a < -1 or a > 1
   */
inline var atanh(const var& x) {
  return make_callback_var(atanh(x.val()), [x](const auto& vi) mutable {
    x.adj() += vi.adj() / (1.0 - x.val() * x.val());
  });
}

/**
 * The inverse hyperbolic tangent function for variables (C99).
 *
 * @tparam Varmat a `var_value` with inner Eigen type
 * @param x The variable.
 * @return Inverse hyperbolic tangent of the variable.
 * @throw std::domain_error if a < -1 or a > 1
 */
template <typename VarMat, require_var_matrix_t<VarMat>* = nullptr>
inline auto atanh(const VarMat& x) {
  return make_callback_var(
      x.val().unaryExpr([](const auto x) { return atanh(x); }),
      [x](const auto& vi) mutable {
        x.adj().array() += vi.adj().array() / (1.0 - x.val().array().square());
      });
}

/**
 * Return the hyperbolic arc tangent of the complex argument.
 *
 * @param[in] z argument
 * @return hyperbolic arc tangent of the argument
 */
inline std::complex<var> atanh(const std::complex<var>& z) {
  return stan::math::internal::complex_atanh(z);
}

}  // namespace math
}  // namespace stan
#endif
