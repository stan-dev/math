#ifndef STAN_MATH_PRIM_FUN_LOG1P_EXP_HPP
#define STAN_MATH_PRIM_FUN_LOG1P_EXP_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Calculates the log of 1 plus the exponential of the specified
 * value without overflow.
 *
 * This function is related to other special functions by:
 *
 * <code>log1p_exp(x) </code>
 *
 * <code> = log1p(exp(a))</code>
 *
 * <code> = log(1 + exp(x))</code>
 *
 * <code> = log_sum_exp(0, x)</code>.
 *
   \f[
   \mbox{log1p\_exp}(x) =
   \begin{cases}
     \ln(1+\exp(x)) & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log1p\_exp}(x)}{\partial x} =
   \begin{cases}
     \frac{\exp(x)}{1+\exp(x)} & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 */
inline double log1p_exp(double a) {
  using std::exp;
  // like log_sum_exp below with b=0.0; prevents underflow
  if (a > 0.0) {
    return a + log1p(exp(-a));
  }
  return log1p(exp(a));
}

/**
 * Structure to wrap log1p_exp() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Natural log of (1 + exp(x)).
 */
struct log1p_exp_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return log1p_exp(x);
  }
};

/**
 * Vectorized version of log1p_exp().
 *
 * @tparam T type of container
 * @param x container
 * @return Natural log of (1 + exp()) applied to each value in x.
 */
template <typename T,
          require_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
          require_not_var_matrix_t<T>* = nullptr>
inline auto log1p_exp(const T& x) {
  return apply_scalar_unary<log1p_exp_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
