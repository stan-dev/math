#ifndef STAN_MATH_PRIM_SCAL_FUN_BETA_HPP
#define STAN_MATH_PRIM_SCAL_FUN_BETA_HPP

#include <stan/math/prim/meta.hpp>
#include <boost/math/tools/promotion.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>

namespace stan {
namespace math {

/**
 * Return the beta function applied to the specified
 * arguments.
 *
 * The beta function is defined for \f$a > 0\f$ and \f$b > 0\f$ by
 *
 * \f$\mbox{B}(a, b) = \frac{\Gamma(a) \Gamma(b)}{\Gamma(a+b)}\f$.
 *
   \f[
   \mbox{beta}(\alpha, \beta) =
   \begin{cases}
     \int_0^1 u^{\alpha - 1} (1 - u)^{\beta - 1} \, du & \mbox{if } \alpha,
 \beta>0 \\[6pt] \textrm{NaN} & \mbox{if } \alpha = \textrm{NaN or } \beta =
 \textrm{NaN} \end{cases} \f]

   \f[
   \frac{\partial\, \mbox{beta}(\alpha, \beta)}{\partial \alpha} =
   \begin{cases}
     \left(\psi(\alpha)-\psi(\alpha+\beta)\right)*\mbox{beta}(\alpha, \beta)
      & \mbox{if } \alpha, \beta>0 \\[6pt]
     \textrm{NaN} & \mbox{if } \alpha = \textrm{NaN or } \beta = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{beta}(\alpha, \beta)}{\partial \beta} =
   \begin{cases}
     \left(\psi(\beta)-\psi(\alpha+\beta)\right)*\mbox{beta}(\alpha, \beta)
      & \mbox{if } \alpha, \beta>0 \\[6pt]
     \textrm{NaN} & \mbox{if } \alpha = \textrm{NaN or } \beta = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam T1 Type of first value.
 * @tparam T2 Type of second value.
 * @param a First value
 * @param b Second value
 * @return Beta function applied to the two values.
 */
template <typename T1, typename T2>
inline return_type_t<T1, T2> beta(const T1 a, const T2 b) {
  using std::exp;
  return exp(lgamma(a) + lgamma(b) - lgamma(a + b));
}

}  // namespace math
}  // namespace stan

#endif
