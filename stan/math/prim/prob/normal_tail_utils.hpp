#ifndef STAN_MATH_PRIM_PROB_NORMAL_TAIL_UTILS_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_TAIL_UTILS_HPP

#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/prob/std_normal_lcdf.hpp>
#include <cmath>

namespace stan {
namespace math {
namespace internal {

template <typename T>
struct normal_tail_terms {
  T log_cdf;
  T mills_ratio;
};

// DLMF 7.12.1 asymptotic expansion of erfcx(x), used above the point where
// forming exp(x^2) directly becomes unsafe.
template <typename T>
inline T erfcx_positive(const T& x) {
  using std::erfc;
  using std::exp;
  if (value_of_rec(x) < 20.0) {
    return exp(x * x) * erfc(x);
  }
  const T inv_two_x_sq = 0.5 / (x * x);
  T term = 1.0;
  T sum = 1.0;
  for (int n = 1; n <= 12; ++n) {
    term *= -(2.0 * n - 1.0) * inv_two_x_sq;
    sum += term;
  }
  return sum * INV_SQRT_PI / x;
}

/**
 * Return the standard normal log CDF and inverse Mills ratio.
 *
 * @tparam T Arithmetic or autodiff scalar type.
 * @param x Standard normal variate.
 * @return `log(Phi(x))` and `phi(x) / Phi(x)`.
 */
template <typename T>
inline normal_tail_terms<T> std_normal_lcdf_and_mills(const T& x) {
  using std::exp;
  using std::log;
  const T erfc_arg = -x * INV_SQRT_TWO;
  if (value_of_rec(erfc_arg) >= 5.0) {
    const T erfcx = erfcx_positive(erfc_arg);
    return {LOG_HALF - erfc_arg * erfc_arg + log(erfcx),
            SQRT_TWO_OVER_SQRT_PI / erfcx};
  }
  const T log_cdf = std_normal_lcdf(x);
  return {log_cdf, exp(NEG_LOG_SQRT_TWO_PI - 0.5 * x * x - log_cdf)};
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
