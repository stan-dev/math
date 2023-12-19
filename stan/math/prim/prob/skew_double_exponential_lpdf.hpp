#ifndef STAN_MATH_PRIM_PROB_SKEW_DOUBLE_EXPONENTIAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_SKEW_DOUBLE_EXPONENTIAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the skew double exponential log probability density function. Given
 * containers of matching sizes, returns the log sum of densities.
 *
 * @tparam T_y type of real parameter.
 * @tparam T_loc type of location parameter.
 * @tparam T_scale type of scale parameter.
 * @tparam T_skewness type of skewness parameter.
 * @param y real parameter
 * @param mu location parameter
 * @param sigma scale parameter
 * @param tau skewness parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if mu is infinite or sigma is nonpositive or tau is
 *  not bound between 0.0 and 1.0
 * @throw std::invalid_argument if container sizes mismatch
 */
template <bool propto, typename T_y, typename T_loc, typename T_scale,
          typename T_skewness,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale, T_skewness>* = nullptr>
return_type_t<T_y, T_loc, T_scale, T_skewness> skew_double_exponential_lpdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma,
    const T_skewness& tau) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_skewness>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  using T_tau_ref = ref_type_if_not_constant_t<T_skewness>;
  static const char* function = "skew_double_exponential_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Shape parameter", sigma, "Skewness parameter",
                         tau);

  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;
  T_tau_ref tau_ref = tau;

  if (size_zero(y, mu, sigma, tau)) {
    return 0.0;
  }
  if (!include_summand<propto, T_y, T_loc, T_scale, T_skewness>::value) {
    return 0.0;
  }

  auto ops_partials
      = make_partials_propagator(y_ref, mu_ref, sigma_ref, tau_ref);

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));
  decltype(auto) tau_val = to_ref(as_value_column_array_or_scalar(tau_ref));

  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive_finite(function, "Scale parameter", sigma_val);
  check_bounded(function, "Skewness parameter", tau_val, 0.0, 1.0);

  const auto& inv_sigma
      = to_ref_if<!is_constant_all<T_scale>::value>(inv(sigma_val));
  const auto& y_m_mu
      = to_ref_if<!is_constant_all<T_y, T_loc>::value>(y_val - mu_val);
  const auto& diff_sign = sign(y_m_mu);

  const auto& diff_sign_smaller_0 = step(-diff_sign);
  const auto& abs_diff_y_mu = fabs(y_m_mu);
  const auto& abs_diff_y_mu_over_sigma = abs_diff_y_mu * inv_sigma;
  const auto& expo = to_ref_if<!is_constant_all<T_skewness>::value>(
      (diff_sign_smaller_0 + diff_sign * tau_val) * abs_diff_y_mu_over_sigma);

  size_t N = max_size(y, mu, sigma, tau);
  T_partials_return logp = -2.0 * sum(expo);

  if (include_summand<propto>::value) {
    logp += N * LOG_TWO;
  }
  if (include_summand<propto, T_scale>::value) {
    logp -= sum(log(sigma_val)) * N / math::size(sigma);
  }
  if (include_summand<propto, T_skewness>::value) {
    logp += sum(log(tau_val) + log1m(tau_val)) * N / math::size(tau);
  }

  if (!is_constant_all<T_y, T_loc>::value) {
    const auto& deriv = to_ref_if<(!is_constant_all<T_y>::value
                                   && !is_constant_all<T_loc>::value)>(
        2.0 * (diff_sign_smaller_0 + diff_sign * tau_val) * diff_sign
        * inv_sigma);
    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials) = -deriv;
    }
    if (!is_constant_all<T_loc>::value) {
      partials<1>(ops_partials) = deriv;
    }
  }
  if (!is_constant_all<T_scale>::value) {
    partials<2>(ops_partials) = -inv_sigma + 2.0 * expo * inv_sigma;
  }
  if (!is_constant_all<T_skewness>::value) {
    edge<3>(ops_partials).partials_
        = inv(tau_val) - inv(1.0 - tau_val)
          + (-1.0 * diff_sign) * 2.0 * abs_diff_y_mu_over_sigma;
  }

  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale, typename T_skewness>
return_type_t<T_y, T_loc, T_scale, T_skewness> skew_double_exponential_lpdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma,
    const T_skewness& tau) {
  return skew_double_exponential_lpdf<false>(y, mu, sigma, tau);
}

}  // namespace math
}  // namespace stan
#endif
