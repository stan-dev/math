#ifndef STAN_MATH_PRIM_PROB_SKEW_NORMAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_SKEW_NORMAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

template <bool propto, typename T_y, typename T_loc, typename T_scale,
          typename T_shape>
return_type_t<T_y, T_loc, T_scale, T_shape> skew_normal_lpdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma, const T_shape& alpha) {
  static const char* function = "skew_normal_lpdf";
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_shape>;

  using std::exp;
  using std::log;

  if (size_zero(y, mu, sigma, alpha)) {
    return 0.0;
  }

  T_partials_return logp(0.0);

  check_not_nan(function, "Random variable", y);
  check_finite(function, "Location parameter", mu);
  check_finite(function, "Shape parameter", alpha);
  check_positive(function, "Scale parameter", sigma);
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma, "Shape paramter", alpha);

  if (!include_summand<propto, T_y, T_loc, T_scale, T_shape>::value) {
    return 0.0;
  }

  operands_and_partials<T_y, T_loc, T_scale, T_shape> ops_partials(y, mu, sigma,
                                                                   alpha);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_loc> mu_vec(mu);
  scalar_seq_view<T_scale> sigma_vec(sigma);
  scalar_seq_view<T_shape> alpha_vec(alpha);
  size_t N = max_size(y, mu, sigma, alpha);

  VectorBuilder<true, T_partials_return, T_scale> inv_sigma(size(sigma));
  VectorBuilder<include_summand<propto, T_scale>::value, T_partials_return,
                T_scale>
      log_sigma(size(sigma));
  for (size_t i = 0; i < size(sigma); i++) {
    inv_sigma[i] = 1.0 / value_of(sigma_vec[i]);
    if (include_summand<propto, T_scale>::value) {
      log_sigma[i] = log(value_of(sigma_vec[i]));
    }
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return mu_dbl = value_of(mu_vec[n]);
    const T_partials_return sigma_dbl = value_of(sigma_vec[n]);
    const T_partials_return alpha_dbl = value_of(alpha_vec[n]);

    const T_partials_return y_minus_mu_over_sigma
        = (y_dbl - mu_dbl) * inv_sigma[n];

    if (include_summand<propto>::value) {
      logp -= HALF_LOG_TWO_PI;
    }
    if (include_summand<propto, T_scale>::value) {
      logp -= log(sigma_dbl);
    }
    if (include_summand<propto, T_y, T_loc, T_scale>::value) {
      logp -= y_minus_mu_over_sigma * y_minus_mu_over_sigma / 2.0;
    }
    logp += log(erfc(-alpha_dbl * y_minus_mu_over_sigma / SQRT_TWO));

    T_partials_return deriv_logerf
        = TWO_OVER_SQRT_PI
          * exp(-alpha_dbl * y_minus_mu_over_sigma / SQRT_TWO * alpha_dbl
                * y_minus_mu_over_sigma / SQRT_TWO)
          / (1 + erf(alpha_dbl * y_minus_mu_over_sigma / SQRT_TWO));
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n]
          += -y_minus_mu_over_sigma / sigma_dbl
             + deriv_logerf * alpha_dbl / (sigma_dbl * SQRT_TWO);
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n]
          += y_minus_mu_over_sigma / sigma_dbl
             + deriv_logerf * -alpha_dbl / (sigma_dbl * SQRT_TWO);
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          += -1.0 / sigma_dbl
             + y_minus_mu_over_sigma * y_minus_mu_over_sigma / sigma_dbl
             - deriv_logerf * y_minus_mu_over_sigma * alpha_dbl
                   / (sigma_dbl * SQRT_TWO);
    }
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge4_.partials_[n]
          += deriv_logerf * y_minus_mu_over_sigma / SQRT_TWO;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
inline return_type_t<T_y, T_loc, T_scale, T_shape> skew_normal_lpdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma, const T_shape& alpha) {
  return skew_normal_lpdf<false>(y, mu, sigma, alpha);
}

}  // namespace math
}  // namespace stan
#endif
