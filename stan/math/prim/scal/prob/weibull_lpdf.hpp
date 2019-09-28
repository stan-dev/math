#ifndef STAN_MATH_PRIM_SCAL_PROB_WEIBULL_LPDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_WEIBULL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the Weibull log probability density for the given
 * location and scale. Given containers of matching sizes, returns the
 * log sum of probability densities.
 *
 * @tparam T_y type of real parameter
 * @tparam T_shape type of shape parameter
 * @tparam T_scale type of scale paramater
 * @param y real parameter
 * @param alpha shape parameter
 * @param sigma scale parameter
 * @return log probability density or log sum of probability densities
 * @throw std::domain_error if y is negative, alpha sigma is nonpositive
 */
template <bool propto, typename T_y, typename T_shape, typename T_scale>
return_type_t<T_y, T_shape, T_scale> weibull_lpdf(const T_y& y,
                                                  const T_shape& alpha,
                                                  const T_scale& sigma) {
  static const char* function = "weibull_lpdf";
  using T_partials_return = partials_return_t<T_y, T_shape, T_scale>;

  using std::log;
  check_finite(function, "Random variable", y);
  check_positive_finite(function, "Shape parameter", alpha);
  check_positive_finite(function, "Scale parameter", sigma);
  check_consistent_sizes(function, "Random variable", y, "Shape parameter",
                         alpha, "Scale parameter", sigma);
  if (size_zero(y, alpha, sigma)) {
    return 0;
  }
  if (!include_summand<propto, T_y, T_shape, T_scale>::value) {
    return 0;
  }

  T_partials_return logp(0);
  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_shape> alpha_vec(alpha);
  scalar_seq_view<T_scale> sigma_vec(sigma);
  size_t N = max_size(y, alpha, sigma);

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    if (y_dbl < 0) {
      return LOG_ZERO;
    }
  }

  VectorBuilder<include_summand<propto, T_shape>::value, T_partials_return,
                T_shape>
      log_alpha(length(alpha));
  for (size_t i = 0; i < length(alpha); i++) {
    if (include_summand<propto, T_shape>::value) {
      log_alpha[i] = log(value_of(alpha_vec[i]));
    }
  }

  VectorBuilder<include_summand<propto, T_y, T_shape>::value, T_partials_return,
                T_y>
      log_y(length(y));
  for (size_t i = 0; i < length(y); i++) {
    if (include_summand<propto, T_y, T_shape>::value) {
      log_y[i] = log(value_of(y_vec[i]));
    }
  }

  VectorBuilder<include_summand<propto, T_shape, T_scale>::value,
                T_partials_return, T_scale>
      log_sigma(length(sigma));
  for (size_t i = 0; i < length(sigma); i++) {
    if (include_summand<propto, T_shape, T_scale>::value) {
      log_sigma[i] = log(value_of(sigma_vec[i]));
    }
  }

  VectorBuilder<include_summand<propto, T_y, T_shape, T_scale>::value,
                T_partials_return, T_scale>
      inv_sigma(length(sigma));
  for (size_t i = 0; i < length(sigma); i++) {
    if (include_summand<propto, T_y, T_shape, T_scale>::value) {
      inv_sigma[i] = 1.0 / value_of(sigma_vec[i]);
    }
  }

  VectorBuilder<include_summand<propto, T_y, T_shape, T_scale>::value,
                T_partials_return, T_y, T_shape, T_scale>
      y_div_sigma_pow_alpha(N);
  for (size_t i = 0; i < N; i++) {
    if (include_summand<propto, T_y, T_shape, T_scale>::value) {
      const T_partials_return y_dbl = value_of(y_vec[i]);
      const T_partials_return alpha_dbl = value_of(alpha_vec[i]);
      y_div_sigma_pow_alpha[i] = pow(y_dbl * inv_sigma[i], alpha_dbl);
    }
  }

  operands_and_partials<T_y, T_shape, T_scale> ops_partials(y, alpha, sigma);
  for (size_t n = 0; n < N; n++) {
    const T_partials_return alpha_dbl = value_of(alpha_vec[n]);
    if (include_summand<propto, T_shape>::value) {
      logp += log_alpha[n];
    }
    if (include_summand<propto, T_y, T_shape>::value) {
      logp += (alpha_dbl - 1.0) * log_y[n];
    }
    if (include_summand<propto, T_shape, T_scale>::value) {
      logp -= alpha_dbl * log_sigma[n];
    }
    if (include_summand<propto, T_y, T_shape, T_scale>::value) {
      logp -= y_div_sigma_pow_alpha[n];
    }

    if (!is_constant_all<T_y>::value) {
      const T_partials_return inv_y = 1.0 / value_of(y_vec[n]);
      ops_partials.edge1_.partials_[n]
          += (alpha_dbl - 1.0) * inv_y
             - alpha_dbl * y_div_sigma_pow_alpha[n] * inv_y;
    }
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge2_.partials_[n]
          += 1.0 / alpha_dbl
             + (1.0 - y_div_sigma_pow_alpha[n]) * (log_y[n] - log_sigma[n]);
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          += alpha_dbl * inv_sigma[n] * (y_div_sigma_pow_alpha[n] - 1.0);
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_shape, typename T_scale>
inline return_type_t<T_y, T_shape, T_scale> weibull_lpdf(const T_y& y,
                                                         const T_shape& alpha,
                                                         const T_scale& sigma) {
  return weibull_lpdf<false>(y, alpha, sigma);
}

}  // namespace math
}  // namespace stan
#endif
