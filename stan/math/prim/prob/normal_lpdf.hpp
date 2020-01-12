#ifndef STAN_MATH_PRIM_PROB_NORMAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of the normal density for the specified scalar(s) given
 * the specified mean(s) and deviation(s). y, mu, or sigma can
 * each be either a scalar or a vector. Any vector inputs
 * must be the same length.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each observation/mean/deviation triple.
 * @tparam T_y Underlying type of scalar in sequence.
 * @tparam T_loc Type of location parameter.
 * @tparam T_scale Type of scale parameter.
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) location parameter(s)
 * for the normal distribution.
 * @param sigma (Sequence of) scale parameters for the normal
 * distribution.
 * @return The log of the product of the densities.
 * @throw std::domain_error if the scale is not positive.
 */
template <bool propto, typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> normal_lpdf(const T_y& y,
                                                      const T_loc& mu,
                                                      const T_scale& sigma) {
  static const char* function = "normal_lpdf";
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;

  using std::log;

  if (size_zero(y, mu, sigma)) {
    return 0.0;
  }

  T_partials_return logp(0.0);

  check_not_nan(function, "Random variable", y);
  check_finite(function, "Location parameter", mu);
  check_positive(function, "Scale parameter", sigma);
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
  if (!include_summand<propto, T_y, T_loc, T_scale>::value) {
    return 0.0;
  }

  operands_and_partials<T_y, T_loc, T_scale> ops_partials(y, mu, sigma);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_loc> mu_vec(mu);
  scalar_seq_view<T_scale> sigma_vec(sigma);
  size_t N = max_size(y, mu, sigma);

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

    const T_partials_return y_minus_mu_over_sigma
        = (y_dbl - mu_dbl) * inv_sigma[n];
    const T_partials_return y_minus_mu_over_sigma_squared
        = y_minus_mu_over_sigma * y_minus_mu_over_sigma;

    static double NEGATIVE_HALF = -0.5;

    if (include_summand<propto>::value) {
      logp += NEG_LOG_SQRT_TWO_PI;
    }
    if (include_summand<propto, T_scale>::value) {
      logp -= log_sigma[n];
    }
    logp += NEGATIVE_HALF * y_minus_mu_over_sigma_squared;

    T_partials_return scaled_diff = inv_sigma[n] * y_minus_mu_over_sigma;
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] -= scaled_diff;
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n] += scaled_diff;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          += -inv_sigma[n] + inv_sigma[n] * y_minus_mu_over_sigma_squared;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> normal_lpdf(const T_y& y,
                                                      const T_loc& mu,
                                                      const T_scale& sigma) {
  return normal_lpdf<false>(y, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
