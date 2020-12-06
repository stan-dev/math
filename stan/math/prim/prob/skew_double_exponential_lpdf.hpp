#ifndef STAN_MATH_PRIM_PROB_SKEW_DOUBLE_EXPONENTIAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_SKEW_DOUBLE_EXPONENTIAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/sign.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
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
 * @param y real parameter
 * @param mu location parameter
 * @param sigma scale parameter
 * @param tau skewness parameter
 * @return log probability density or log sum of probability densities
 * @throw std::domain_error if y is nan, mu is infinite, or sigma is nonpositive
 * @throw std::invalid_argument if container sizes mismatch
 */
template <bool propto, typename T_y, typename T_loc, typename T_scale, typename T_skewness>
return_type_t<T_y, T_loc, T_scale, T_skewness> skew_double_exponential_lpdf(
  const T_y& y, const T_loc& mu, const T_scale& sigma, const T_skewness& tau) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_skewness>;
  using T_partials_array = Eigen::Array<T_partials_return, Eigen::Dynamic, 1>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_mu_ref = ref_type_if_t<!is_constant<T_loc>::value, T_loc>;
  using T_sigma_ref = ref_type_if_t<!is_constant<T_scale>::value, T_scale>;
  using T_tau_ref = ref_type_if_t<!is_constant<T_skewness>::value, T_skewness>;
  static const char* function = "skew_double_exponential_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Shape parameter", sigma, "Skewness parameter", tau);

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

  T_partials_return logp(0.0);
  operands_and_partials<T_y_ref, T_mu_ref, T_sigma_ref, T_tau_ref> ops_partials(
      y_ref, mu_ref, sigma_ref, tau_ref);

  const auto& y_col = as_column_vector_or_scalar(y_ref);
  const auto& mu_col = as_column_vector_or_scalar(mu_ref);
  const auto& sigma_col = as_column_vector_or_scalar(sigma_ref);
  const auto& tau_col = as_column_vector_or_scalar(tau_ref);

  const auto& y_arr = as_array_or_scalar(y_col);
  const auto& mu_arr = as_array_or_scalar(mu_col);
  const auto& sigma_arr = as_array_or_scalar(sigma_col);
  const auto& tau_arr = as_column_vector_or_scalar(tau_col);

  ref_type_t<decltype(value_of(y_arr))> y_val = value_of(y_arr);
  ref_type_t<decltype(value_of(mu_arr))> mu_val = value_of(mu_arr);
  ref_type_t<decltype(value_of(sigma_arr))> sigma_val = value_of(sigma_arr);
  ref_type_t<decltype(value_of(sigma_arr))> tau_val = value_of(tau_arr);

  check_finite(function, "Random variable", y);
  check_finite(function, "Location parameter", mu);
  check_positive_finite(function, "Scale parameter", sigma);
  check_positive_finite(function, "Skewness parameter", tau);

  size_t N = max_size(y, mu, sigma);

  logp += log(tau) + log1m(tau) - log(sigma);
  if (y < mu) {
    logp -= 2 * (1 - tau) * (mu - y) / sigma;
  }
  else {
    logp -= 2 * tau * (y - mu)  / sigma;
  }


  //  real skew_double_exponential_lpdf(real y, real mu, real sigma, real tau) {
  //    return log(tau) + log1m(tau)
  //           - log(sigma)
  //           - 2 * ((y < mu) ? (1 - tau) * (mu - y) : tau * (y - mu)) / sigma;
  //  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale, typename T_skewness>
return_type_t<T_y, T_loc, T_scale, T_skewness> skew_double_exponential_lpdf(
  const T_y& y, const T_loc& mu, const T_scale& sigma, const T_skewness& tau) {
  return skew_double_exponential_lpdf<false>(y, mu, sigma, tau);
}

}  // namespace math
}  // namespace stan
#endif
