#ifndef STAN_MATH_PRIM_SCAL_PROB_NORMAL_LPDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_NORMAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <cmath>
#include <utility>

namespace stan {
namespace math {

/**
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
auto normal_lpdf(T_y&& y, T_loc&& mu, T_scale&& sigma) {
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
  const auto sigma_size = length(sigma);
  const size_t N = max_size(y, mu, sigma);
  const scalar_seq_view<T_y> y_vec(std::forward<T_y>(y));
  const scalar_seq_view<T_loc> mu_vec(std::forward<T_loc>(mu));
  const scalar_seq_view<T_scale> sigma_vec(std::forward<T_scale>(sigma));

  VectorBuilder<true, T_partials_return, T_scale> inv_sigma(sigma_size);
  VectorBuilder<include_summand<propto, T_scale>::value, T_partials_return,
                T_scale>
      log_sigma(sigma_size);
  for (size_t i = 0; i < sigma_size; i++) {
    if (include_summand<propto, T_scale>::value) {
      log_sigma[i] = log(value_of(sigma_vec[i]));
    } else {
      inv_sigma[i] = 1.0 / value_of(sigma_vec[i]);
    }
  }

  for (size_t n = 0; n < N; n++) {
    const auto y_dbl = value_of(y_vec[n]);
    const auto mu_dbl = value_of(mu_vec[n]);

    const auto y_scaled = (y_dbl - mu_dbl) * inv_sigma[n];
    const auto y_scaled_sq = y_scaled * y_scaled;

    static double NEGATIVE_HALF = -0.5;

    if (include_summand<propto>::value) {
      logp += NEG_LOG_SQRT_TWO_PI;
    }
    if (include_summand<propto, T_scale>::value) {
      logp -= log_sigma[n];
    }
    if (include_summand<propto, T_y, T_loc, T_scale>::value) {
      logp += NEGATIVE_HALF * y_scaled_sq;
    }

    T_partials_return scaled_diff = inv_sigma[n] * y_scaled;
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] -= scaled_diff;
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n] += scaled_diff;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          += -inv_sigma[n] + inv_sigma[n] * y_scaled_sq;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale>
inline auto normal_lpdf(T_y&& y, T_loc&& mu, T_scale&& sigma) {
  return normal_lpdf<false>(std::forward<T_y>(y), std::forward<T_loc>(mu),
                            std::forward<T_scale>(sigma));
}

}  // namespace math
}  // namespace stan
#endif
