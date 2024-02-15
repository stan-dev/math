#ifndef STAN_MATH_PRIM_PROB_SKEW_DOUBLE_EXPONENTIAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_SKEW_DOUBLE_EXPONENTIAL_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the skew double exponential log complementary cumulative density
 * function. Given containers of matching sizes, returns the log sum of
 * probabilities.
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
template <typename T_y, typename T_loc, typename T_scale, typename T_skewness,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale, T_skewness>* = nullptr>
return_type_t<T_y, T_loc, T_scale, T_skewness> skew_double_exponential_lccdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma,
    const T_skewness& tau) {
  static constexpr const char* function = "skew_double_exponential_lccdf";
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_skewness>;
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Shape parameter", sigma, "Skewness parameter",
                         tau);
  auto&& y_ref = to_ref(y);
  auto&& mu_ref = to_ref(mu);
  auto&& sigma_ref = to_ref(sigma);
  auto&& tau_ref = to_ref(tau);

  auto&& y_val = as_value_array_or_scalar(y_ref);
  auto&& mu_val = as_value_array_or_scalar(mu_ref);
  auto&& sigma_val = as_value_array_or_scalar(sigma_ref);
  auto&& tau_val = as_value_array_or_scalar(tau_ref);

  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive_finite(function, "Scale parameter", sigma_val);
  check_bounded(function, "Skewness parameter", tau_val, 0.0, 1.0);
  if (size_zero(y, mu, sigma, tau)) {
    return 0.0;
  }

  auto ops_partials
      = make_partials_propagator(y_ref, mu_ref, sigma_ref, tau_ref);

  scalar_seq_view<std::decay_t<decltype(y_val)>> y_vec(y_val);
  scalar_seq_view<std::decay_t<decltype(mu_val)>> mu_vec(mu_val);
  scalar_seq_view<std::decay_t<decltype(sigma_val)>> sigma_vec(sigma_val);
  scalar_seq_view<std::decay_t<decltype(tau_val)>> tau_vec(tau_val);

  const auto N = max_size(y, mu, sigma, tau);
  auto inv_sigma_val = to_ref(inv(sigma_val));
  scalar_seq_view<decltype(inv_sigma_val)> inv_sigma(inv_sigma_val);

  T_partials_return cdf_log(0.0);
  for (int i = 0; i < N; ++i) {
    const T_partials_return y_dbl = y_vec[i];
    const T_partials_return mu_dbl = mu_vec[i];
    const T_partials_return sigma_dbl = sigma_vec[i];
    const T_partials_return tau_dbl = tau_vec[i];

    const T_partials_return y_m_mu = y_dbl - mu_dbl;
    const T_partials_return diff_sign = sign(y_m_mu);
    const T_partials_return diff_sign_smaller_0 = step(-diff_sign);
    const T_partials_return abs_diff_y_mu = fabs(y_m_mu);
    const T_partials_return abs_diff_y_mu_over_sigma
        = abs_diff_y_mu * inv_sigma[i];
    const T_partials_return expo = (diff_sign_smaller_0 + diff_sign * tau_dbl)
                                   * abs_diff_y_mu_over_sigma;
    const T_partials_return inv_exp_2_expo_tau
        = inv(exp(2.0 * (tau_dbl - 1.0) * y_m_mu * inv_sigma[i]) - tau_dbl);
    const T_partials_return tau_m1_tau = (tau_dbl - 1.0) * tau_dbl;

    const T_partials_return rep_deriv
        = y_dbl < mu_dbl ? 2.0 * tau_m1_tau * inv_sigma[i] * inv_exp_2_expo_tau
                         : -2.0 * inv_sigma[i] * tau_dbl;
    const T_partials_return sig_deriv
        = y_dbl < mu_dbl ? -2.0 * tau_m1_tau * y_m_mu * inv_sigma[i]
                               * inv_sigma[i] * inv_exp_2_expo_tau
                         : 2.0 * inv_sigma[i] * expo;
    const T_partials_return skew_deriv
        = y_dbl < mu_dbl ? -(sigma_dbl + 2.0 * tau_dbl * y_m_mu * diff_sign)
                               * inv_sigma[i] * inv_exp_2_expo_tau
                         : 1.0 / (tau_dbl - 1.0) - 2.0 * inv_sigma[i] * y_m_mu;

    if (y_dbl <= mu_dbl) {
      cdf_log += log1m(tau_dbl * exp(-2.0 * expo));
    } else {
      cdf_log += log1m(tau_dbl) - 2.0 * expo;
    }

    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials)[i] += rep_deriv;
    }
    if (!is_constant_all<T_loc>::value) {
      partials<1>(ops_partials)[i] -= rep_deriv;
    }
    if (!is_constant_all<T_scale>::value) {
      partials<2>(ops_partials)[i] += sig_deriv;
    }
    if (!is_constant_all<T_skewness>::value) {
      partials<3>(ops_partials)[i] += skew_deriv;
    }
  }
  return ops_partials.build(cdf_log);
}
}  // namespace math
}  // namespace stan
#endif
