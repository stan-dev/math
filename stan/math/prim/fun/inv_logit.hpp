#ifndef STAN_MATH_PRIM_FUN_INV_LOGIT_HPP
#define STAN_MATH_PRIM_FUN_INV_LOGIT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

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
    if (a < LOG_EPSILON) {
      return exp_a;
    }
    return exp_a / (1 + exp_a);
  }
  return inv(1 + exp(-a));
}

/**
 * Structure to wrap inv_logit() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
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
 *
 * @tparam T type of container
 * @param x container
 * @return Inverse logit applied to each value in x.
 */
template <
    typename T, require_not_var_matrix_t<T>* = nullptr,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto inv_logit(const T& x) {
  return apply_scalar_unary<inv_logit_fun, T>::apply(x);
}

// TODO(Tadej): Eigen is introducing their implementation logistic() of this
// in 3.4. Use that once we switch to Eigen 3.4

}  // namespace math
}  // namespace stan

#endif
