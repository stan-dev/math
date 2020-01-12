#ifndef STAN_MATH_PRIM_PROB_SKEW_NORMAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_SKEW_NORMAL_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/owens_t.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
return_type_t<T_y, T_loc, T_scale, T_shape> skew_normal_lccdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma, const T_shape& alpha) {
  static const char* function = "skew_normal_lccdf";
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_shape>;

  T_partials_return ccdf_log(0.0);

  if (size_zero(y, mu, sigma, alpha)) {
    return ccdf_log;
  }

  check_not_nan(function, "Random variable", y);
  check_finite(function, "Location parameter", mu);
  check_not_nan(function, "Scale parameter", sigma);
  check_positive(function, "Scale parameter", sigma);
  check_finite(function, "Shape parameter", alpha);
  check_not_nan(function, "Shape parameter", alpha);
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma, "Shape paramter", alpha);

  operands_and_partials<T_y, T_loc, T_scale, T_shape> ops_partials(y, mu, sigma,
                                                                   alpha);

  using std::exp;
  using std::log;

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_loc> mu_vec(mu);
  scalar_seq_view<T_scale> sigma_vec(sigma);
  scalar_seq_view<T_shape> alpha_vec(alpha);
  size_t N = max_size(y, mu, sigma, alpha);

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return mu_dbl = value_of(mu_vec[n]);
    const T_partials_return sigma_dbl = value_of(sigma_vec[n]);
    const T_partials_return alpha_dbl = value_of(alpha_vec[n]);
    const T_partials_return alpha_dbl_sq = alpha_dbl * alpha_dbl;
    const T_partials_return diff = (y_dbl - mu_dbl) / sigma_dbl;
    const T_partials_return diff_sq = diff * diff;
    const T_partials_return scaled_diff = diff / SQRT_TWO;
    const T_partials_return scaled_diff_sq = diff_sq * 0.5;
    const T_partials_return ccdf_log_
        = 1.0 - 0.5 * erfc(-scaled_diff) + 2 * owens_t(diff, alpha_dbl);

    ccdf_log += log(ccdf_log_);

    const T_partials_return deriv_erfc
        = INV_SQRT_TWO_PI * exp(-scaled_diff_sq) / sigma_dbl;
    const T_partials_return deriv_owens
        = erf(alpha_dbl * scaled_diff) * exp(-scaled_diff_sq)
          / SQRT_TWO_OVER_SQRT_PI / (-TWO_PI) / sigma_dbl;
    const T_partials_return rep_deriv
        = (-2.0 * deriv_owens + deriv_erfc) / ccdf_log_;

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] -= rep_deriv;
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n] += rep_deriv;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n] += rep_deriv * diff;
    }
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge4_.partials_[n]
          -= -2.0 * exp(-0.5 * diff_sq * (1.0 + alpha_dbl_sq))
             / ((1 + alpha_dbl_sq) * TWO_PI) / ccdf_log_;
    }
  }
  return ops_partials.build(ccdf_log);
}

}  // namespace math
}  // namespace stan
#endif
