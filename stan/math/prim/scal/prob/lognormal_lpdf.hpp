#ifndef STAN_MATH_PRIM_SCAL_PROB_LOGNORMAL_LPDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_LOGNORMAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
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

// LogNormal(y|mu, sigma)  [y >= 0;  sigma > 0]
template <bool propto, typename T_y, typename T_loc, typename T_scale>
inline auto lognormal_lpdf(const T_y& y, const T_loc& mu,
                           const T_scale& sigma) {
  using T_partials = partials_return_t<T_y, T_loc, T_scale>;
  T_partials logp(0);
  using T_return = return_type_t<T_y, T_loc, T_scale>;

  using std::log;
  static const char* function = "lognormal_lpdf";
  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_finite(function, "Location parameter", mu);
  check_positive_finite(function, "Scale parameter", sigma);
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);

  const scalar_seq_view<T_y> y_vec(y);
  const scalar_seq_view<T_loc> mu_vec(mu);
  const scalar_seq_view<T_scale> sigma_vec(sigma);
  const size_t N = max_size(y, mu, sigma);
  operands_and_partials<T_y, T_loc, T_scale> ops_partials(y, mu, sigma);
  if (size_zero(y, mu, sigma)) {
    return ops_partials.build(logp);
  }

  for (size_t n = 0; n < length(y); n++) {
    if (value_of(y_vec[n]) <= 0) {
      return ops_partials.build(T_partials(LOG_ZERO));
    }
  }

  if (include_summand<propto>::value) {
    logp += N * NEG_LOG_SQRT_TWO_PI;
  }

  T_partials log_y = 0;
  T_partials logy_m_mu(0);
  T_partials inv_sigma_sq = 0;
  T_partials inv_sigma = 0;
  T_partials log_sigma = 0;
  for (size_t n = 0; n < N; n++) {
    const T_partials mu_dbl = value_of(mu_vec[n]);
    const auto y_val = value_of(y_vec[n]);
    const auto sigma_val = value_of(sigma_vec[n]);
    if (include_summand<propto, T_y, T_loc, T_scale>::value) {
      log_y = log(y_val);
      inv_sigma = 1 / sigma_val;
      inv_sigma_sq = inv_sigma * inv_sigma;
      logy_m_mu = log_y - mu_dbl;
    }

    T_partials logy_m_mu_sq = logy_m_mu * logy_m_mu;
    T_partials logy_m_mu_div_sigma(0);
    if (!is_constant_all<T_y, T_loc, T_scale>::value) {
      logy_m_mu_div_sigma = logy_m_mu * inv_sigma_sq;
    }

    if (include_summand<propto, T_scale>::value) {
      log_sigma = log(sigma_val);
      logp -= log_sigma;
    }
    if (include_summand<propto, T_y>::value) {
      logp -= log_y;
    }
    if (include_summand<propto, T_y, T_loc, T_scale>::value) {
      logp -= 0.5 * logy_m_mu_sq * inv_sigma_sq;
    }

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n]
          -= (1 + logy_m_mu_div_sigma) * (1 / y_val);
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n] += logy_m_mu_div_sigma;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          += (logy_m_mu_div_sigma * logy_m_mu - 1) * inv_sigma;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale>
inline auto lognormal_lpdf(const T_y& y, const T_loc& mu,
                           const T_scale& sigma) {
  return lognormal_lpdf<false>(y, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
