#ifndef STAN_MATH_PRIM_FUN_LOG1P_EXP_HPP
#define STAN_MATH_PRIM_FUN_LOG1P_EXP_HPP

#include <stanh/prim/fun/log1p.hpp>
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
  if (a > 0.0)
    return a + log1p(exp(-a));
  return log1p(exp(a));
}

}  // namespace math
}  // namespace stan

#endif
#ifndef STAN_MATH_PRIM_FUN_LOG1P_EXP_HPP
#define STAN_MATH_PRIM_FUN_LOG1P_EXP_HPP

#include <stanh/prim/vectorize/apply_scalar_unary.hpp>
#include <stanh/prim/fun/log1p_exp.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap log1m_exp() so that it can be vectorized.
 * @param x Variable.
 * @tparam T Variable type.
 * @return Natural log of (1 + exp(x)).
 */
struct log1p_exp_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return log1p_exp(x);
  }
};

/**
 * Vectorized version of log1m_exp().
 * @param x Container.
 * @tparam T Container type.
 * @return Natural log of (1 + exp()) applied to each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<log1p_exp_fun, T>::return_t log1p_exp(
    const T& x) {
  return apply_scalar_unary<log1p_exp_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
#ifndef STAN_MATH_PRIM_FUN_LOG1P_EXP_HPP
#define STAN_MATH_PRIM_FUN_LOG1P_EXP_HPP

#include <stanh/prim/fun/log1p.hpp>
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
  if (a > 0.0)
    return a + log1p(exp(-a));
  return log1p(exp(a));
}

}  // namespace math
}  // namespace stan

#endif
