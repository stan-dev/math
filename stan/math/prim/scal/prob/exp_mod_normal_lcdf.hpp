#ifndef STAN_MATH_PRIM_SCAL_PROB_EXP_MOD_NORMAL_LCDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_EXP_MOD_NORMAL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/is_inf.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale, typename T_inv_scale>
inline auto exp_mod_normal_lcdf(T_y&& y, T_loc&& mu,
                                T_scale&& sigma,
                                T_inv_scale&& lambda) {
  using T_partials = partials_return_t<T_y, T_loc, T_scale, T_inv_scale>;
  T_partials cdf_log(0.0);
  using std::exp;
  using std::log;

  static const char* function = "exp_mod_normal_lcdf";
  check_not_nan(function, "Random variable", y);
  check_finite(function, "Location parameter", mu);
  check_not_nan(function, "Scale parameter", sigma);
  check_positive_finite(function, "Scale parameter", sigma);
  check_positive_finite(function, "Inv_scale parameter", lambda);
  check_not_nan(function, "Inv_scale parameter", lambda);
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma, "Inv_scale paramter",
                         lambda);

  const double sqrt_pi = std::sqrt(pi());
  const scalar_seq_view<T_y> y_vec(y);
  const scalar_seq_view<T_loc> mu_vec(mu);
  const scalar_seq_view<T_scale> sigma_vec(sigma);
  const scalar_seq_view<T_inv_scale> lambda_vec(lambda);
  const size_t N = max_size(y, mu, sigma, lambda);
  operands_and_partials<T_y, T_loc, T_scale, T_inv_scale> ops_partials(
      y, mu, sigma, lambda);
  if (size_zero(y, mu, sigma, lambda)) {
    return ops_partials.build(cdf_log);
  }
  for (size_t n = 0; n < N; n++) {
    if (is_inf(y_vec[n])) {
      if (y_vec[n] < 0.0) {
        return ops_partials.build(T_partials(negative_infinity()));
      } else {
        return ops_partials.build(cdf_log);
      }
    }

    const auto y_dbl = value_of(y_vec[n]);
    const auto mu_dbl = value_of(mu_vec[n]);
    const auto sigma_dbl = value_of(sigma_vec[n]);
    const auto lambda_dbl = value_of(lambda_vec[n]);
    const auto u = lambda_dbl * (y_dbl - mu_dbl);
    const auto v = lambda_dbl * sigma_dbl;
    const auto v_sq = v * v;
    const auto scaled_diff = (y_dbl - mu_dbl) / (SQRT_2 * sigma_dbl);
    const auto scaled_diff_sq = scaled_diff * scaled_diff;
    const auto erf_calc1 = 0.5 * (1 + erf(u / (v * SQRT_2)));
    const auto erf_calc2 = 0.5 * (1 + erf(u / (v * SQRT_2) - v / SQRT_2));
    const auto deriv_1 = lambda_dbl * exp(0.5 * v_sq - u) * erf_calc2;
    const auto deriv_2
        = SQRT_2 / sqrt_pi * 0.5
          * exp(0.5 * v_sq
                - (-scaled_diff + (v / SQRT_2)) * (-scaled_diff + (v / SQRT_2))
                - u)
          / sigma_dbl;
    const auto deriv_3
        = SQRT_2 / sqrt_pi * 0.5 * exp(-scaled_diff_sq) / sigma_dbl;

    const auto denom = erf_calc1 - erf_calc2 * exp(0.5 * v_sq - u);
    const auto cdf_ = erf_calc1 - exp(-u + v_sq * 0.5) * (erf_calc2);

    cdf_log += log(cdf_);

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] += (deriv_1 - deriv_2 + deriv_3) / denom;
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n]
          += (-deriv_1 + deriv_2 - deriv_3) / denom;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          += (-deriv_1 * v - deriv_3 * scaled_diff * SQRT_2
              - deriv_2 * sigma_dbl * SQRT_2
                    * (-SQRT_2 * 0.5
                           * (-lambda_dbl + scaled_diff * SQRT_2 / sigma_dbl)
                       - SQRT_2 * lambda_dbl))
             / denom;
    }
    if (!is_constant_all<T_inv_scale>::value) {
      ops_partials.edge4_.partials_[n]
          += exp(0.5 * v_sq - u)
             * (SQRT_2 / sqrt_pi * 0.5 * sigma_dbl
                    * exp(-(v / SQRT_2 - scaled_diff)
                          * (v / SQRT_2 - scaled_diff))
                - (v * sigma_dbl + mu_dbl - y_dbl) * erf_calc2)
             / denom;
    }
  }
  return ops_partials.build(cdf_log);
}

}  // namespace math
}  // namespace stan
#endif
