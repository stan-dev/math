#ifndef STAN_MATH_REV_CORE_OPERATOR_UNARY_PLUS_HPP
#define STAN_MATH_REV_CORE_OPERATOR_UNARY_PLUS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/callback_vari.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_nan.hpp>

namespace stan {
namespace math {

/**
 * Unary plus operator for variables (C++).
 *
 * The function simply returns its input, because
 *
 * \f$\frac{d}{dx} +x = \frac{d}{dx} x = 1\f$.
 *
 * The effect of unary plus on a built-in C++ scalar type is
 * integer promotion.  Because variables are all
 * double-precision floating point already, promotion is
 * not necessary.
 *
 \f[
 \mbox{operator+}(x) =
 \begin{cases}
 x & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
 \textrm{NaN} & \mbox{if } x = \textrm{NaN}
 \end{cases}
 \f]

 \f[
 \frac{\partial\, \mbox{operator+}(x)}{\partial x} =
 \begin{cases}
 1 & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
 \textrm{NaN} & \mbox{if } x = \textrm{NaN}
 \end{cases}
 \f]
 *
 * @param a Argument variable.
 * @return The input reference.
 */
inline var operator+(const var& a) { return a; }

}  // namespace math
}  // namespace stan
#endif
