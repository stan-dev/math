#ifndef STAN_MATH_REV_FUN_HYPOT_HPP
#define STAN_MATH_REV_FUN_HYPOT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/hypot.hpp>

namespace stan {
namespace math {

/**
 * Returns the length of the hypotenuse of a right triangle
 * with sides of the specified lengths (C99).
 *
 * The partial derivatives are given by
 *
 * \f$\frac{\partial}{\partial x} \sqrt{x^2 + y^2} = \frac{x}{\sqrt{x^2 +
 * y^2}}\f$, and
 *
 * \f$\frac{\partial}{\partial y} \sqrt{x^2 + y^2} = \frac{y}{\sqrt{x^2 +
 * y^2}}\f$.
 *
 * @param[in] a Length of first side.
 * @param[in] b Length of second side.
 * @return Length of hypotenuse.
 */
inline var hypot(const var& a, const var& b) {
  return make_callback_var(hypot(a.val(), b.val()), [a, b](auto& vi) mutable {
    a.adj() += vi.adj() * a.val() / vi.val();
    b.adj() += vi.adj() * b.val() / vi.val();
  });
}

/**
 * Returns the length of the hypotenuse of a right triangle
 * with sides of the specified lengths (C99).
 *
 * The derivative is
 *
 * \f$\frac{d}{d x} \sqrt{x^2 + c^2} = \frac{x}{\sqrt{x^2 + c^2}}\f$.
 *
 * @param[in] a Length of first side.
 * @param[in] b Length of second side.
 * @return Length of hypotenuse.
 */
inline var hypot(const var& a, double b) {
  return make_callback_var(hypot(a.val(), b), [a](auto& vi) mutable {
    a.adj() += vi.adj() * a.val() / vi.val();
  });
}

/**
 * Returns the length of the hypotenuse of a right triangle
 * with sides of the specified lengths (C99).
 *
 * The derivative is
 *
 * \f$\frac{d}{d y} \sqrt{c^2 + y^2} = \frac{y}{\sqrt{c^2 + y^2}}\f$.
 *
   \f[
   \mbox{hypot}(x, y) =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < 0 \text{ or } y < 0 \\
     \sqrt{x^2+y^2} & \mbox{if } x, y\geq 0 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{hypot}(x, y)}{\partial x} =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < 0 \text{ or } y < 0 \\
     \frac{x}{\sqrt{x^2+y^2}} & \mbox{if } x, y\geq 0 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{hypot}(x, y)}{\partial y} =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < 0 \text{ or } y < 0 \\
     \frac{y}{\sqrt{x^2+y^2}} & \mbox{if } x, y\geq 0 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param[in] a Length of first side.
 * @param[in] b Length of second side.
 * @return Length of hypotenuse.
 */
inline var hypot(double a, const var& b) {
  return make_callback_var(hypot(b.val(), a), [b](auto& vi) mutable {
    b.adj() += vi.adj() * b.val() / vi.val();
  });
}

}  // namespace math
}  // namespace stan
#endif
