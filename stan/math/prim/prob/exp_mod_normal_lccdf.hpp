#ifndef STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/is_inf.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale, typename T_inv_scale>
return_type_t<T_y, T_loc, T_scale, T_inv_scale> exp_mod_normal_lccdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma,
    const T_inv_scale& lambda) {
  static const char* function = "exp_mod_normal_lccdf";
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_inv_scale>;

  T_partials_return ccdf_log(0.0);
  if (size_zero(y, mu, sigma, lambda)) {
    return ccdf_log;
  }

  check_not_nan(function, "Random variable", y);
  check_finite(function, "Location parameter", mu);
  check_not_nan(function, "Scale parameter", sigma);
  check_positive_finite(function, "Scale parameter", sigma);
  check_positive_finite(function, "Inv_scale parameter", lambda);
  check_not_nan(function, "Inv_scale parameter", lambda);
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma, "Inv_scale paramter",
                         lambda);

  operands_and_partials<T_y, T_loc, T_scale, T_inv_scale> ops_partials(
      y, mu, sigma, lambda);

  using std::exp;
  using std::log;

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_loc> mu_vec(mu);
  scalar_seq_view<T_scale> sigma_vec(sigma);
  scalar_seq_view<T_inv_scale> lambda_vec(lambda);
  size_t N = max_size(y, mu, sigma, lambda);

  for (size_t n = 0, size_y = size(y); n < size_y; n++) {
    if (is_inf(y_vec[n])) {
      return ops_partials.build(y_vec[n] > 0 ? negative_infinity() : 0);
    }
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return mu_dbl = value_of(mu_vec[n]);
    const T_partials_return sigma_dbl = value_of(sigma_vec[n]);
    const T_partials_return lambda_dbl = value_of(lambda_vec[n]);
    const T_partials_return inv_sigma = inv(sigma_dbl);
    const T_partials_return diff = y_dbl - mu_dbl;
    const T_partials_return u = lambda_dbl * diff;
    const T_partials_return v = lambda_dbl * sigma_dbl;
    const T_partials_return v_over_sqrt_two = v * INV_SQRT_TWO;
    const T_partials_return scaled_diff = diff * INV_SQRT_TWO * inv_sigma;
    const T_partials_return erf_calc
        = 0.5 * (1 + erf(u / (v * SQRT_TWO) - v_over_sqrt_two));
    const T_partials_return exp_term = exp(0.5 * square(v) - u);

    const T_partials_return deriv_1 = lambda_dbl * exp_term * erf_calc;
    const T_partials_return deriv_2
        = INV_SQRT_TWO_PI * exp_term
          * exp(-square(-scaled_diff + v_over_sqrt_two)) * inv_sigma;
    const T_partials_return deriv_3
        = INV_SQRT_TWO_PI * exp(-square(scaled_diff)) * inv_sigma;

    const T_partials_return ccdf_n
        = 0.5 - 0.5 * erf(u / (v * SQRT_TWO)) + exp_term * erf_calc;

    ccdf_log += log(ccdf_n);

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n]
          -= (deriv_1 - deriv_2 + deriv_3) / ccdf_n;
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n]
          += (deriv_1 - deriv_2 + deriv_3) / ccdf_n;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          += (deriv_1 * v + deriv_3 * scaled_diff * SQRT_TWO
              - deriv_2 * (v + SQRT_TWO * scaled_diff))
             / ccdf_n;
    }
    if (!is_constant_all<T_inv_scale>::value) {
      ops_partials.edge4_.partials_[n]
          -= exp_term
             * (INV_SQRT_TWO_PI * sigma_dbl
                    * exp(-square(v_over_sqrt_two - scaled_diff))
                - (v * sigma_dbl - diff) * erf_calc)
             / ccdf_n;
    }
  }

  return ops_partials.build(ccdf_log);
}

}  // namespace math
}  // namespace stan
#endif
