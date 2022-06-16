#ifndef STAN_MATH_PRIM_FUN_LOG1M_EXP_HPP
#define STAN_MATH_PRIM_FUN_LOG1M_EXP_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/expm1.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

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
  if (a > 0) {
    return NOT_A_NUMBER;
  } else if (a > -0.693147) {
    return log(-expm1(a));  // 0.693147 ~= log(2)
  } else {
    return log1m(exp(a));
  }
}

/**
 * Structure to wrap log1m_exp() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Natural log of (1 - exp(x)).
 */
struct log1m_exp_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    return log1m_exp(x);
  }
};

/**
 * Vectorized version of log1m_exp().
 *
 * @tparam T type of container
 * @param x container
 * @return Natural log of (1 - exp()) applied to each value in x.
 */
template <
    typename T, require_not_var_matrix_t<T>* = nullptr,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto log1m_exp(const T& x) {
  return apply_scalar_unary<log1m_exp_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
