#ifndef STAN_MATH_PRIM_SCAL_PROB_LOGNORMAL_LCDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_LOGNORMAL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale>
inline auto lognormal_lcdf(T_y&& y, T_loc&& mu, T_scale&& sigma) {
  using T_partials = partials_return_t<T_y, T_loc, T_scale>;
  T_partials cdf_log = 0.0;

  using std::exp;
  using std::log;

  static const char* function = "lognormal_lcdf";
  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_finite(function, "Location parameter", mu);
  check_positive_finite(function, "Scale parameter", sigma);
  const scalar_seq_view<T_y> y_vec(y);
  const scalar_seq_view<T_loc> mu_vec(mu);
  const scalar_seq_view<T_scale> sigma_vec(sigma);
  const size_t N = max_size(y, mu, sigma);
  const double sqrt_pi = std::sqrt(pi());
  operands_and_partials<T_y, T_loc, T_scale> ops_partials(y, mu, sigma);
  if (size_zero(y, mu, sigma)) {
    return ops_partials.build(cdf_log);
  }

  for (size_t i = 0; i < stan::length(y); i++) {
    if (value_of(y_vec[i]) == 0.0) {
      return ops_partials.build(T_partials(negative_infinity()));
    }
  }

  const double log_half = std::log(0.5);

  for (size_t n = 0; n < N; n++) {
    const T_partials y_dbl = value_of(y_vec[n]);
    const T_partials mu_dbl = value_of(mu_vec[n]);
    const T_partials sigma_dbl = value_of(sigma_vec[n]);
    const T_partials scaled_diff = (log(y_dbl) - mu_dbl) / (sigma_dbl * SQRT_2);
    const T_partials rep_deriv
        = SQRT_2 / sqrt_pi * exp(-scaled_diff * scaled_diff) / sigma_dbl;

    const T_partials erfc_calc = erfc(-scaled_diff);
    cdf_log += log_half + log(erfc_calc);

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] += rep_deriv / erfc_calc / y_dbl;
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n] -= rep_deriv / erfc_calc;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          -= rep_deriv * scaled_diff * SQRT_2 / erfc_calc;
    }
  }
  return ops_partials.build(cdf_log);
}

}  // namespace math
}  // namespace stan
#endif
