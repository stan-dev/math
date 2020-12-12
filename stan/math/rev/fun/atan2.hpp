#ifndef STAN_MATH_REV_FUN_ATAN2_HPP
#define STAN_MATH_REV_FUN_ATAN2_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the principal value of the arc tangent, in radians, of
 * the first variable divided by the second (cmath).
 *
 * The partial derivatives are defined by
 *
 * \f$ \frac{\partial}{\partial x} \arctan \frac{x}{y} = \frac{y}{x^2 + y^2}\f$,
 * and
 *
 * \f$ \frac{\partial}{\partial y} \arctan \frac{x}{y} = \frac{-x}{x^2 +
 * y^2}\f$.
 *
 * @param a Numerator variable.
 * @param b Denominator variable.
 * @return The arc tangent of the fraction, in radians.
 */
inline var atan2(const var& a, const var& b) {
  return make_callback_var(
      std::atan2(a.val(), b.val()), [a, b](const auto& vi) mutable {
        double a_sq_plus_b_sq = (a.val() * a.val()) + (b.val() * b.val());
        a.adj() += vi.adj_ * b.val() / a_sq_plus_b_sq;
        b.adj() -= vi.adj_ * a.val() / a_sq_plus_b_sq;
      });
}

/**
 * Return the principal value of the arc tangent, in radians, of
 * the first variable divided by the second scalar (cmath).
 *
 * The derivative with respect to the variable is
 *
 * \f$ \frac{d}{d x} \arctan \frac{x}{c} = \frac{c}{x^2 + c^2}\f$.
 *
 * @param a Numerator variable.
 * @param b Denominator scalar.
 * @return The arc tangent of the fraction, in radians.
 */
inline var atan2(const var& a, double b) {
  return make_callback_var(
      std::atan2(a.val(), b), [a, b](const auto& vi) mutable {
        double a_sq_plus_b_sq = (a.val() * a.val()) + (b * b);
        a.adj() += vi.adj_ * b / a_sq_plus_b_sq;
      });
}

/**
 * Return the principal value of the arc tangent, in radians, of
 * the first scalar divided by the second variable (cmath).
 *
 * The derivative with respect to the variable is
 *
 * \f$ \frac{\partial}{\partial y} \arctan \frac{c}{y} = \frac{-c}{c^2 +
 y^2}\f$.
 *
 *
   \f[
   \mbox{atan2}(x, y) =
   \begin{cases}
     \arctan\left(\frac{x}{y}\right) & \mbox{if } -\infty\leq x \leq \infty,
 -\infty\leq y \leq \infty \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or
 } y = \textrm{NaN} \end{cases} \f]

   \f[
   \frac{\partial\, \mbox{atan2}(x, y)}{\partial x} =
   \begin{cases}
     \frac{y}{x^2+y^2} & \mbox{if } -\infty\leq x\leq \infty, -\infty\leq y \leq
 \infty \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{atan2}(x, y)}{\partial y} =
   \begin{cases}
     -\frac{x}{x^2+y^2} & \mbox{if } -\infty\leq x\leq \infty, -\infty\leq y
 \leq \infty \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y =
 \textrm{NaN} \end{cases} \f]
 *
 * @param a Numerator scalar.
 * @param b Denominator variable.
 * @return The arc tangent of the fraction, in radians.
 */
inline var atan2(double a, const var& b) {
  return make_callback_var(
      std::atan2(a, b.val()), [a, b](const auto& vi) mutable {
        double a_sq_plus_b_sq = (a * a) + (b.val() * b.val());
        b.adj() -= vi.adj_ * a / a_sq_plus_b_sq;
      });
}

}  // namespace math
}  // namespace stan
#endif
