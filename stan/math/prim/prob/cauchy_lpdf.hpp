#ifndef STAN_MATH_PRIM_PROB_CAUCHY_LPDF_HPP
#define STAN_MATH_PRIM_PROB_CAUCHY_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of the Cauchy density for the specified scalar(s) given
 * the specified location parameter(s) and scale parameter(s). y,
 * mu, or sigma can each either be scalar a vector.  Any vector
 * inputs must be the same length.
 *
 * <p> The result log probability is defined to be the sum of
 * the log probabilities for each observation/mu/sigma triple.
 *
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) location(s).
 * @param sigma (Sequence of) scale(s).
 * @return The log of the product of densities.
 * @tparam T_y Type of scalar outcome.
 * @tparam T_loc Type of location.
 * @tparam T_scale Type of scale.
 */
template <bool propto, typename T_y, typename T_loc, typename T_scale>
return_type_t<T_y, T_loc, T_scale> cauchy_lpdf(const T_y& y, const T_loc& mu,
                                               const T_scale& sigma) {
  static const char* function = "cauchy_lpdf";
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;

  if (size_zero(y, mu, sigma)) {
    return 0.0;
  }

  T_partials_return logp(0.0);

  check_not_nan(function, "Random variable", y);
  check_finite(function, "Location parameter", mu);
  check_positive_finite(function, "Scale parameter", sigma);
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);

  if (!include_summand<propto, T_y, T_loc, T_scale>::value) {
    return 0.0;
  }

  using std::log;

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_loc> mu_vec(mu);
  scalar_seq_view<T_scale> sigma_vec(sigma);
  size_t N = max_size(y, mu, sigma);

  VectorBuilder<true, T_partials_return, T_scale> inv_sigma(size(sigma));
  VectorBuilder<true, T_partials_return, T_scale> sigma_squared(size(sigma));
  VectorBuilder<include_summand<propto, T_scale>::value, T_partials_return,
                T_scale>
      log_sigma(size(sigma));
  for (size_t i = 0; i < size(sigma); i++) {
    const T_partials_return sigma_dbl = value_of(sigma_vec[i]);
    inv_sigma[i] = 1.0 / sigma_dbl;
    sigma_squared[i] = sigma_dbl * sigma_dbl;
    if (include_summand<propto, T_scale>::value) {
      log_sigma[i] = log(sigma_dbl);
    }
  }

  operands_and_partials<T_y, T_loc, T_scale> ops_partials(y, mu, sigma);

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return mu_dbl = value_of(mu_vec[n]);

    const T_partials_return y_minus_mu = y_dbl - mu_dbl;
    const T_partials_return y_minus_mu_squared = y_minus_mu * y_minus_mu;
    const T_partials_return y_minus_mu_over_sigma = y_minus_mu * inv_sigma[n];
    const T_partials_return y_minus_mu_over_sigma_squared
        = y_minus_mu_over_sigma * y_minus_mu_over_sigma;

    if (include_summand<propto>::value) {
      logp -= LOG_PI;
    }
    if (include_summand<propto, T_scale>::value) {
      logp -= log_sigma[n];
    }
    logp -= log1p(y_minus_mu_over_sigma_squared);

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n]
          -= 2 * y_minus_mu / (sigma_squared[n] + y_minus_mu_squared);
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n]
          += 2 * y_minus_mu / (sigma_squared[n] + y_minus_mu_squared);
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          += (y_minus_mu_squared - sigma_squared[n]) * inv_sigma[n]
             / (sigma_squared[n] + y_minus_mu_squared);
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> cauchy_lpdf(const T_y& y,
                                                      const T_loc& mu,
                                                      const T_scale& sigma) {
  return cauchy_lpdf<false>(y, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
