#ifndef STAN_MATH_PRIM_SCAL_PROB_CAUCHY_LPDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_CAUCHY_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/log1p.hpp>
#include <cmath>
#include <utility>

namespace stan {
namespace math {

/**
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
inline auto cauchy_lpdf(T_y&& y, T_loc&& mu, T_scale&& sigma) {
  static const char* function = "cauchy_lpdf";
  using T_partials = partials_return_t<T_y, T_loc, T_scale>;
  T_partials logp(0.0);

  using std::log;

  check_not_nan(function, "Random variable", y);
  check_finite(function, "Location parameter", mu);
  check_positive_finite(function, "Scale parameter", sigma);
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);

  const scalar_seq_view<T_y> y_vec(y);
  const scalar_seq_view<T_loc> mu_vec(mu);
  const scalar_seq_view<T_scale> sigma_vec(sigma);
  const size_t N = max_size(y, mu, sigma);
  operands_and_partials<T_y, T_loc, T_scale> ops_partials(y, mu, sigma);
  if (!include_summand<propto, T_y, T_loc, T_scale>::value) {
    return ops_partials.build(logp);
  } else if (size_zero(y, mu, sigma)) {
    return ops_partials.build(logp);
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials y_dbl = value_of(y_vec[n]);
    const T_partials mu_dbl = value_of(mu_vec[n]);
    const T_partials sigma_dbl = value_of(sigma_vec[n]);
    const T_partials inv_sigma = 1.0 / sigma_dbl;
    const T_partials sigma_squared = sigma_dbl * sigma_dbl;
    const T_partials y_minus_mu = y_dbl - mu_dbl;
    const T_partials y_minus_mu_squared = y_minus_mu * y_minus_mu;
    const T_partials y_minus_mu_over_sigma = y_minus_mu * inv_sigma;
    const T_partials y_minus_mu_over_sigma_squared
        = y_minus_mu_over_sigma * y_minus_mu_over_sigma;

    if (include_summand<propto>::value) {
      logp += NEG_LOG_PI;
    }
    if (include_summand<propto, T_scale>::value) {
      logp -= log(sigma_dbl);
    }
    if (include_summand<propto, T_y, T_loc, T_scale>::value) {
      logp -= log1p(y_minus_mu_over_sigma_squared);
    }

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n]
          -= 2 * y_minus_mu / (sigma_squared + y_minus_mu_squared);
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n]
          += 2 * y_minus_mu / (sigma_squared + y_minus_mu_squared);
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          += (y_minus_mu_squared - sigma_squared) * inv_sigma
             / (sigma_squared + y_minus_mu_squared);
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale>
inline auto cauchy_lpdf(T_y&& y, T_loc&& mu, T_scale&& sigma) {
  return cauchy_lpdf<false>(std::forward<T_y>(y), std::forward<T_loc>(mu), std::forward<T_scale>(sigma));
}

}  // namespace math
}  // namespace stan
#endif
