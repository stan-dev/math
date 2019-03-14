#ifndef STAN_MATH_PRIM_FUN_INV_LOGIT_HPP
#define STAN_MATH_PRIM_FUN_INV_LOGIT_HPP

#include <stan/math/prim/fun/constants.hpp>
#include <cmath>
#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>






namespace stan {
namespace math {

/**
 * Returns the inverse logit function applied to the argument.
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
 * @param a Argument.
 * @return Inverse logit of argument.
 */
inline double inv_logit(double a) {
  using std::exp;
  if (a < 0) {
    double exp_a = exp(a);
    if (a < LOG_EPSILON)
      return exp_a;
    return exp_a / (1 + exp_a);
  }
  return 1 / (1 + exp(-a));
}













/**
 * Structure to wrap inv_logit() so that it can be vectorized.
 * @param x Variable.
 * @tparam T Variable type.
 * @return Inverse logit of x.
 */
struct inv_logit_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return inv_logit(x);
  }
};

/**
 * Vectorized version of inv_logit().
 * @param x Container.
 * @tparam T Container type.
 * @return Inverse logit applied to each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<inv_logit_fun, T>::return_t inv_logit(
    const T& x) {
  return apply_scalar_unary<inv_logit_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
