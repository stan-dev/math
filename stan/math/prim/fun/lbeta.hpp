#ifndef STAN_MATH_PRIM_SCAL_FUN_LBETA_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LBETA_HPP

#include <limits>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/lgamma_stirling.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/lgamma_stirling_diff.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/err/check_nonnegative.hpp>
#include <limits>

namespace stan {
namespace math {

/**
 * Return the log of the beta function applied to the specified
 * arguments.
 *
 * The beta function is defined for \f$a > 0\f$ and \f$b > 0\f$ by
 *
 * \f$\mbox{B}(a, b) = \frac{\Gamma(a) \Gamma(b)}{\Gamma(a+b)}\f$.
 *
 * This function returns its log,
 *
 * \f$\log \mbox{B}(a, b) = \log \Gamma(a) + \log \Gamma(b) - \log
 \Gamma(a+b)\f$.
 *
 * See stan::math::lgamma() for the double-based and stan::math for the
 * variable-based log Gamma function.
 * This function is numerically more stable than naive evaluation via lgamma
 *
   \f[
   \mbox{lbeta}(\alpha, \beta) =
   \begin{cases}
     \ln\int_0^1 u^{\alpha - 1} (1 - u)^{\beta - 1} \, du & \mbox{if } \alpha,
 \beta>0 \\[6pt] \textrm{NaN} & \mbox{if } \alpha = \textrm{NaN or } \beta =
 \textrm{NaN} \end{cases} \f]

   \f[
   \frac{\partial\, \mbox{lbeta}(\alpha, \beta)}{\partial \alpha} =
   \begin{cases}
     \Psi(\alpha)-\Psi(\alpha+\beta) & \mbox{if } \alpha, \beta>0 \\[6pt]
     \textrm{NaN} & \mbox{if } \alpha = \textrm{NaN or } \beta = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{lbeta}(\alpha, \beta)}{\partial \beta} =
   \begin{cases}
     \Psi(\beta)-\Psi(\alpha+\beta) & \mbox{if } \alpha, \beta>0 \\[6pt]
     \textrm{NaN} & \mbox{if } \alpha = \textrm{NaN or } \beta = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a First value
 * @param b Second value
 * @return Log of the beta function applied to the two values.
 * @tparam T1 Type of first value.
 * @tparam T2 Type of second value.
 */
template <typename T1, typename T2>
return_type_t<T1, T2> lbeta(const T1 a, const T2 b) {
  typedef return_type_t<T1, T2> T_ret;

  if (is_nan(value_of_rec(a)) || is_nan(value_of_rec(b))) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  static const char* function = "lbeta";
  check_nonnegative(function, "first argument", a);
  check_nonnegative(function, "second argument", b);
  T_ret x;  // x is the smaller of the two
  T_ret y;
  if (a < b) {
    x = a;
    y = b;
  } else {
    x = b;
    y = a;
  }

  // Special cases
  if (x == 0) {
    return std::numeric_limits<double>::infinity();
  } else if (is_inf(y)) {
    return -std::numeric_limits<double>::infinity();
  }

  // For large x or y, separate the lgamma values into Stirling approximations
  // and appropriate corrections. The Stirling approximations allow for
  // analytic simplifaction and the corrections are added later.
  //
  // The overall approach is inspired by the code in R, where the algorithm is
  // credited to W. Fullerton of Los Alamos Scientific Laboratory
  if (y < lgamma_stirling_diff_useful) {
    // both small
    return lgamma(x) + lgamma(y) - lgamma(x + y);
  } else if (x < lgamma_stirling_diff_useful) {
    // y large, x small
    T_ret stirling_diff = lgamma_stirling_diff(y) - lgamma_stirling_diff(x + y);
    T_ret log_x_y = log(x + y);
    T_ret stirling = (y - 0.5) * log1p(-x / (x + y)) + x * (1 - log_x_y);
    return stirling + lgamma(x) + stirling_diff;
  } else {
    // both large
    T_ret stirling_diff = lgamma_stirling_diff(x) + lgamma_stirling_diff(y)
                          - lgamma_stirling_diff(x + y);
    T_ret stirling = (x - 0.5) * log(x / (x + y)) + y * log1p(-x / (x + y))
                     + 0.5 * (LOG_TWO_PI - log(y));
    return stirling + stirling_diff;
  }
}

}  // namespace math
}  // namespace stan

#endif
