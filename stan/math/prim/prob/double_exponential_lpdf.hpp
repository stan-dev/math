#ifndef STAN_MATH_PRIM_PROB_DOUBLE_EXPONENTIAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_DOUBLE_EXPONENTIAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/sign.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the double exponential log probability density function. Given
 * containers of matching sizes, returns the log sum of densities.
 *
 * @tparam T_y type of real parameter.
 * @tparam T_loc type of location parameter.
 * @tparam T_scale type of scale parameter.
 * @param y real parameter
 * @param mu location parameter
 * @param sigma scale parameter
 * @return log probability density or log sum of probability densities
 * @throw std::domain_error if y is nan, mu is infinite, or sigma is nonpositive
 * @throw std::invalid_argument if container sizes mismatch
 */
template <bool propto, typename T_y, typename T_loc, typename T_scale>
return_type_t<T_y, T_loc, T_scale> double_exponential_lpdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma) {
  static const char* function = "double_exponential_lpdf";
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;

  using std::fabs;
  using std::log;

  if (size_zero(y, mu, sigma)) {
    return 0.0;
  }

  T_partials_return logp(0.0);
  check_finite(function, "Random variable", y);
  check_finite(function, "Location parameter", mu);
  check_positive_finite(function, "Scale parameter", sigma);
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Shape parameter", sigma);

  if (!include_summand<propto, T_y, T_loc, T_scale>::value) {
    return 0.0;
  }

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_loc> mu_vec(mu);
  scalar_seq_view<T_scale> sigma_vec(sigma);
  size_t N = max_size(y, mu, sigma);
  operands_and_partials<T_y, T_loc, T_scale> ops_partials(y, mu, sigma);

  VectorBuilder<include_summand<propto, T_y, T_loc, T_scale>::value,
                T_partials_return, T_scale>
      inv_sigma(size(sigma));
  VectorBuilder<!is_constant_all<T_scale>::value, T_partials_return, T_scale>
      inv_sigma_squared(size(sigma));
  VectorBuilder<include_summand<propto, T_scale>::value, T_partials_return,
                T_scale>
      log_sigma(size(sigma));
  for (size_t i = 0; i < size(sigma); i++) {
    const T_partials_return sigma_dbl = value_of(sigma_vec[i]);
    inv_sigma[i] = 1.0 / sigma_dbl;
    if (include_summand<propto, T_scale>::value) {
      log_sigma[i] = log(value_of(sigma_vec[i]));
    }
    if (!is_constant_all<T_scale>::value) {
      inv_sigma_squared[i] = inv_sigma[i] * inv_sigma[i];
    }
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return mu_dbl = value_of(mu_vec[n]);

    const T_partials_return y_m_mu = y_dbl - mu_dbl;
    const T_partials_return fabs_y_m_mu = fabs(y_m_mu);

    if (include_summand<propto>::value) {
      logp -= LOG_TWO;
    }
    if (include_summand<propto, T_scale>::value) {
      logp -= log_sigma[n];
    }
    logp -= fabs_y_m_mu * inv_sigma[n];

    T_partials_return sign_y_m_mu_times_inv_sigma(0);
    if (!is_constant_all<T_y, T_loc>::value) {
      sign_y_m_mu_times_inv_sigma = sign(y_m_mu) * inv_sigma[n];
    }
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] -= sign_y_m_mu_times_inv_sigma;
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n] += sign_y_m_mu_times_inv_sigma;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          += -inv_sigma[n] + fabs_y_m_mu * inv_sigma_squared[n];
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale>
return_type_t<T_y, T_loc, T_scale> double_exponential_lpdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma) {
  return double_exponential_lpdf<false>(y, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
