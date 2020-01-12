#ifndef STAN_MATH_PRIM_PROB_LOGNORMAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_LOGNORMAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

// LogNormal(y|mu, sigma)  [y >= 0;  sigma > 0]
template <bool propto, typename T_y, typename T_loc, typename T_scale>
return_type_t<T_y, T_loc, T_scale> lognormal_lpdf(const T_y& y, const T_loc& mu,
                                                  const T_scale& sigma) {
  static const char* function = "lognormal_lpdf";
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;

  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_finite(function, "Location parameter", mu);
  check_positive_finite(function, "Scale parameter", sigma);
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
  if (size_zero(y, mu, sigma)) {
    return 0;
  }

  T_partials_return logp(0);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_loc> mu_vec(mu);
  scalar_seq_view<T_scale> sigma_vec(sigma);
  size_t N = max_size(y, mu, sigma);

  for (size_t n = 0; n < size(y); n++) {
    if (value_of(y_vec[n]) <= 0) {
      return LOG_ZERO;
    }
  }

  operands_and_partials<T_y, T_loc, T_scale> ops_partials(y, mu, sigma);

  using std::log;

  VectorBuilder<include_summand<propto, T_scale>::value, T_partials_return,
                T_scale>
      log_sigma(size(sigma));
  if (include_summand<propto, T_scale>::value) {
    for (size_t n = 0; n < size(sigma); n++) {
      log_sigma[n] = log(value_of(sigma_vec[n]));
    }
  }

  VectorBuilder<include_summand<propto, T_y, T_loc, T_scale>::value,
                T_partials_return, T_scale>
      inv_sigma(size(sigma));
  VectorBuilder<include_summand<propto, T_y, T_loc, T_scale>::value,
                T_partials_return, T_scale>
      inv_sigma_sq(size(sigma));
  if (include_summand<propto, T_y, T_loc, T_scale>::value) {
    for (size_t n = 0; n < size(sigma); n++) {
      inv_sigma[n] = 1 / value_of(sigma_vec[n]);
    }
  }
  if (include_summand<propto, T_y, T_loc, T_scale>::value) {
    for (size_t n = 0; n < size(sigma); n++) {
      inv_sigma_sq[n] = inv_sigma[n] * inv_sigma[n];
    }
  }

  VectorBuilder<include_summand<propto, T_y, T_loc, T_scale>::value,
                T_partials_return, T_y>
      log_y(size(y));
  if (include_summand<propto, T_y, T_loc, T_scale>::value) {
    for (size_t n = 0; n < size(y); n++) {
      log_y[n] = log(value_of(y_vec[n]));
    }
  }

  VectorBuilder<!is_constant_all<T_y>::value, T_partials_return, T_y> inv_y(
      size(y));
  if (!is_constant_all<T_y>::value) {
    for (size_t n = 0; n < size(y); n++) {
      inv_y[n] = 1 / value_of(y_vec[n]);
    }
  }

  if (include_summand<propto>::value) {
    logp += N * NEG_LOG_SQRT_TWO_PI;
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials_return mu_dbl = value_of(mu_vec[n]);

    T_partials_return logy_m_mu(0);
    if (include_summand<propto, T_y, T_loc, T_scale>::value) {
      logy_m_mu = log_y[n] - mu_dbl;
    }

    T_partials_return logy_m_mu_sq = logy_m_mu * logy_m_mu;
    T_partials_return logy_m_mu_div_sigma(0);
    if (!is_constant_all<T_y, T_loc, T_scale>::value) {
      logy_m_mu_div_sigma = logy_m_mu * inv_sigma_sq[n];
    }

    if (include_summand<propto, T_scale>::value) {
      logp -= log_sigma[n];
    }
    if (include_summand<propto, T_y>::value) {
      logp -= log_y[n];
    }
    if (include_summand<propto, T_y, T_loc, T_scale>::value) {
      logp -= 0.5 * logy_m_mu_sq * inv_sigma_sq[n];
    }

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] -= (1 + logy_m_mu_div_sigma) * inv_y[n];
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n] += logy_m_mu_div_sigma;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          += (logy_m_mu_div_sigma * logy_m_mu - 1) * inv_sigma[n];
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> lognormal_lpdf(const T_y& y,
                                                         const T_loc& mu,
                                                         const T_scale& sigma) {
  return lognormal_lpdf<false>(y, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
