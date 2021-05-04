#ifndef STAN_MATH_PRIM_FUN_LGAMMA_STIRLING_DIFF_HPP
#define STAN_MATH_PRIM_FUN_LGAMMA_STIRLING_DIFF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/is_nan.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/lgamma_stirling.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

constexpr double lgamma_stirling_diff_useful = 10;

/**
 * Return the difference between log of the gamma function and its Stirling
 * approximation.
 * This is useful to stably compute log of ratios of gamma functions with large
 * arguments where the Stirling approximation allows for analytic solution
 * and the (small) differences can be added afterwards.
 * This is for example used in the implementation of lbeta.
 *
 * The function will return correct value for all arguments, but using it can
 * lead to a loss of precision when x < lgamma_stirling_diff_useful.
 *
   \f[
   \mbox{lgamma_stirling_diff}(x) =
    \log(\Gamma(x)) - \frac{1}{2} \log(2\pi) +
    (x-\frac{1}{2})*\log(x) - x
   \f]

 *
 * @tparam T type of value
 * @param x value
 * @return Difference between lgamma(x) and its Stirling approximation.
 */
template <typename T>
return_type_t<T> lgamma_stirling_diff(const T x) {
  using T_ret = return_type_t<T>;

  if (is_nan(value_of_rec(x))) {
    return NOT_A_NUMBER;
  }
  check_nonnegative("lgamma_stirling_diff", "argument", x);

  if (x == 0) {
    return INFTY;
  }
  if (value_of(x) < lgamma_stirling_diff_useful) {
    return lgamma(x) - lgamma_stirling(x);
  }

  // Using the Stirling series as expressed in formula 5.11.1. at
  // https://dlmf.nist.gov/5.11
  constexpr double stirling_series[]{
      0.0833333333333333333333333,   -0.00277777777777777777777778,
      0.000793650793650793650793651, -0.000595238095238095238095238,
      0.000841750841750841750841751, -0.00191752691752691752691753,
      0.00641025641025641025641026,  -0.0295506535947712418300654};

  constexpr int n_stirling_terms = 6;
  T_ret result(0.0);
  T_ret multiplier = inv(x);
  T_ret inv_x_squared = square(multiplier);
  for (int n = 0; n < n_stirling_terms; n++) {
    if (n > 0) {
      multiplier *= inv_x_squared;
    }
    result += stirling_series[n] * multiplier;
  }
  return result;
}

}  // namespace math
}  // namespace stan

#endif
