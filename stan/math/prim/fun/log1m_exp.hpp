#ifndef STAN_MATH_PRIM_FUN_LOG1M_EXP_HPP
#define STAN_MATH_PRIM_FUN_LOG1M_EXP_HPP

#include <stanh/prim/fun/expm1.hpp>
#include <stanh/prim/fun/log1m.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

/**
 * Calculates the natural logarithm of one minus the exponential
 * of the specified value without overflow,
 *
 * <p><code>log1m_exp(x) = log(1-exp(x))</code>
 *
 * This function is only defined for x < 0
 *

   \f[
   \mbox{log1m\_exp}(x) =
   \begin{cases}
     \ln(1-\exp(x)) & \mbox{if } x < 0 \\
     \textrm{NaN} & \mbox{if } x \geq 0\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{asinh}(x)}{\partial x} =
   \begin{cases}
     -\frac{\exp(x)}{1-\exp(x)} & \mbox{if } x < 0 \\
     \textrm{NaN} & \mbox{if } x \geq 0\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

 * @param[in] a Argument.
 * @return natural logarithm of one minus the exponential of the
 * argument.
 *
 */
inline double log1m_exp(double a) {
  using std::exp;
  using std::log;
  if (a >= 0)
    return std::numeric_limits<double>::quiet_NaN();
  else if (a > -0.693147)
    return log(-expm1(a));  // 0.693147 ~= log(2)
  else
    return log1m(exp(a));
}

}  // namespace math
}  // namespace stan

#endif
#ifndef STAN_MATH_PRIM_FUN_LOG1M_EXP_HPP
#define STAN_MATH_PRIM_FUN_LOG1M_EXP_HPP

#include <stanh/prim/vectorize/apply_scalar_unary.hpp>
#include <stanh/prim/fun/log1m_exp.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap log1m_exp() so it can be vectorized.
 * @param x Variable.
 * @tparam T Variable type.
 * @return Natural log of (1 - exp(x)).
 */
struct log1m_exp_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return log1m_exp(x);
  }
};

/**
 * Vectorized version of log1m_exp().
 * @param x Container.
 * @tparam T Container type.
 * @return Natural log of (1 - exp()) applied to each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<log1m_exp_fun, T>::return_t log1m_exp(
    const T& x) {
  return apply_scalar_unary<log1m_exp_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
#ifndef STAN_MATH_PRIM_FUN_LOG1M_EXP_HPP
#define STAN_MATH_PRIM_FUN_LOG1M_EXP_HPP

#include <stanh/prim/fun/expm1.hpp>
#include <stanh/prim/fun/log1m.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

/**
 * Calculates the natural logarithm of one minus the exponential
 * of the specified value without overflow,
 *
 * <p><code>log1m_exp(x) = log(1-exp(x))</code>
 *
 * This function is only defined for x < 0
 *

   \f[
   \mbox{log1m\_exp}(x) =
   \begin{cases}
     \ln(1-\exp(x)) & \mbox{if } x < 0 \\
     \textrm{NaN} & \mbox{if } x \geq 0\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{asinh}(x)}{\partial x} =
   \begin{cases}
     -\frac{\exp(x)}{1-\exp(x)} & \mbox{if } x < 0 \\
     \textrm{NaN} & \mbox{if } x \geq 0\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

 * @param[in] a Argument.
 * @return natural logarithm of one minus the exponential of the
 * argument.
 *
 */
inline double log1m_exp(double a) {
  using std::exp;
  using std::log;
  if (a >= 0)
    return std::numeric_limits<double>::quiet_NaN();
  else if (a > -0.693147)
    return log(-expm1(a));  // 0.693147 ~= log(2)
  else
    return log1m(exp(a));
}

}  // namespace math
}  // namespace stan

#endif
