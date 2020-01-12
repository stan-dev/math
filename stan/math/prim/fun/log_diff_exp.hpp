#ifndef STAN_MATH_PRIM_FUN_LOG_DIFF_EXP_HPP
#define STAN_MATH_PRIM_FUN_LOG_DIFF_EXP_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/log1m_exp.hpp>
#include <stan/math/prim/fun/constants.hpp>

namespace stan {
namespace math {

/**
 * The natural logarithm of the difference of the natural exponentiation
 * of x and the natural exponentiation of y
 *
 * This function is only defined for x >= y
 *
 *
   \f[
   \mbox{log\_diff\_exp}(x, y) =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < y\\
     \ln(\exp(x)-\exp(y)) & \mbox{if } x \geq y \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log\_diff\_exp}(x, y)}{\partial x} =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x \leq y\\
     \frac{\exp(x)}{\exp(x)-\exp(y)} & \mbox{if } x > y \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log\_diff\_exp}(x, y)}{\partial y} =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x \leq y\\
     -\frac{\exp(y)}{\exp(x)-\exp(y)} & \mbox{if } x > y \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 */
template <typename T1, typename T2>
inline return_type_t<T1, T2> log_diff_exp(const T1 x, const T2 y) {
  if (x <= y) {
    return (x < INFTY && x == y) ? NEGATIVE_INFTY : NOT_A_NUMBER;
  }
  return x + log1m_exp(y - x);
}

}  // namespace math
}  // namespace stan

#endif
