#ifndef STAN_MATH_OPENCL_PRIM_FUN_INV_SQUARE_HPP
#define STAN_MATH_OPENCL_PRIM_FUN_INV_SQUARE_HPP

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/fun/constants.hpp>

namespace stan {
namespace math {
/**
 * Returns the inverse logit function applied to the kernel generator expression.
 *
 * The inverse logit function is defined by
 *
 * \f$\mbox{logit}^{-1}(x) = \frac{1}{1 + \exp(-x)}\f$.
 *
 * This function can be used to implement the inverse link function
 * for logistic regression.
 *
 * The inverse to this function is <code>logit</code>.
 *
 \f[
 \mbox{inv\_logit}(y) =
 \begin{cases}
 \mbox{logit}^{-1}(y) & \mbox{if } -\infty\leq y \leq \infty \\[6pt]
 \textrm{NaN} & \mbox{if } y = \textrm{NaN}
 \end{cases}
 \f]

 \f[
 \frac{\partial\, \mbox{inv\_logit}(y)}{\partial y} =
 \begin{cases}
 \frac{\partial\, \mbox{logit}^{-1}(y)}{\partial y} & \mbox{if } -\infty\leq
 y\leq \infty \\[6pt] \textrm{NaN} & \mbox{if } y = \textrm{NaN} \end{cases} \f]

 \f[
 \mbox{logit}^{-1}(y) = \frac{1}{1 + \exp(-y)}
 \f]

 \f[
 \frac{\partial \, \mbox{logit}^{-1}(y)}{\partial y} =
 \frac{\exp(y)}{(\exp(y)+1)^2} \f]
 *
 * @param a input kernel generator expression.
 * @return inverse logit of each value in x.
 */
template <typename T_a,
          typename = require_all_valid_expressions_and_none_scalar_t<T_a>>
inline auto inv_logit(T_a&& a) {  // NOLINT
  return select((a < 0.0),
                select((a < LOG_EPSILON),
                  exp(a),
                  elewise_division(exp(a), 1.0 + exp(a))),
                elewise_division(1.0, 1.0 + exp(-a)));
}
}  // namespace math
}  // namespace stan

#endif