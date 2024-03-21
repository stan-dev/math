#ifndef STAN_MATH_PRIM_FUN_LBETA_HPP
#define STAN_MATH_PRIM_FUN_LBETA_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <stan/math/prim/fun/is_inf.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/lgamma_stirling_diff.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>

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
 * This function is numerically more stable than naive evaluation via lgamma.
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
 * @tparam T1 type of first value
 * @tparam T2 type of second value
 * @param a First value
 * @param b Second value
 * @return Log of the beta function applied to the two values.
 */
template <typename T1, typename T2, require_all_arithmetic_t<T1, T2>* = nullptr>
return_type_t<T1, T2> lbeta(const T1 a, const T2 b) {
  using T_ret = return_type_t<T1, T2>;

  if (is_any_nan(a, b)) {
    return NOT_A_NUMBER;
  }

  static constexpr const char* function = "lbeta";
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
    return INFTY;
  }
  if (is_inf(y)) {
    return NEGATIVE_INFTY;
  }

  // For large x or y, separate the lgamma values into Stirling approximations
  // and appropriate corrections. The Stirling approximations allow for
  // analytic simplification and the corrections are added later.
  //
  // The overall approach is inspired by the code in R, where the algorithm is
  // credited to W. Fullerton of Los Alamos Scientific Laboratory
  if (y < lgamma_stirling_diff_useful) {
    // both small
    return lgamma(x) + lgamma(y) - lgamma(x + y);
  }
  T_ret x_over_xy = x / (x + y);
  if (x < lgamma_stirling_diff_useful) {
    // y large, x small
    T_ret stirling_diff = lgamma_stirling_diff(y) - lgamma_stirling_diff(x + y);
    T_ret stirling = (y - 0.5) * log1m(x_over_xy) + x * (1 - log(x + y));
    return stirling + lgamma(x) + stirling_diff;
  }

  // both large
  T_ret stirling_diff = lgamma_stirling_diff(x) + lgamma_stirling_diff(y)
                        - lgamma_stirling_diff(x + y);
  T_ret stirling = (x - 0.5) * log(x_over_xy) + y * log1m(x_over_xy)
                   + HALF_LOG_TWO_PI - 0.5 * log(y);
  return stirling + stirling_diff;
}

/**
 * Enables the vectorized application of the lbeta function,
 * when the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return lbeta function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto lbeta(const T1& a, const T2& b) {
  return apply_scalar_binary(
      a, b, [&](const auto& c, const auto& d) { return lbeta(c, d); });
}

}  // namespace math
}  // namespace stan

#endif
