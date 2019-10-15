#ifndef STAN_MATH_PRIM_SCAL_FUN_LOG_SUM_EXP_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LOG_SUM_EXP_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/log1p_exp.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <boost/math/tools/promotion.hpp>

namespace stan {
namespace math {

/**
 * Calculates the log sum of exponetials without overflow.
 *
 * \f$\log (\exp(a) + \exp(b)) = m + \log(\exp(a-m) + \exp(b-m))\f$,
 *
 * where \f$m = max(a, b)\f$.
 *
 *
   \f[
   \mbox{log\_sum\_exp}(x, y) =
   \begin{cases}
     \ln(\exp(x)+\exp(y)) & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log\_sum\_exp}(x, y)}{\partial x} =
   \begin{cases}
     \frac{\exp(x)}{\exp(x)+\exp(y)} & \mbox{if } -\infty\leq x, y \leq \infty
 \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log\_sum\_exp}(x, y)}{\partial y} =
   \begin{cases}
     \frac{\exp(y)}{\exp(x)+\exp(y)} & \mbox{if } -\infty\leq x, y \leq \infty
 \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a the first variable
 * @param b the second variable
 */
template <typename T1, typename T2>
inline return_type_t<T1, T2> log_sum_exp(const T2& a, const T1& b) {
  if (a == NEGATIVE_INFTY) {
    return b;
  }
  if (a == INFTY && b == INFTY) {
    return INFTY;
  }
  if (a > b) {
    return a + log1p_exp(b - a);
  }
  return b + log1p_exp(a - b);
}

}  // namespace math
}  // namespace stan

#endif
