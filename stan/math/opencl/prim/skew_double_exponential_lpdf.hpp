#ifndef STAN_MATH_OPENCL_PRIM_SKEW_DOUBLE_EXPONENTIAL_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_SKEW_DOUBLE_EXPONENTIAL_LPDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns the log PMF of the skew double exponential distribution. If
 * containers are supplied, returns the log sum of the probabilities.
 *
 * @tparam T_y_cl type of dependent variable
 * @tparam T_loc_cl type of location parameter
 * @tparam T_scale_cl type of scale parameter
 * @tparam T_skewness_cl type of inverse scale parameter
 * @param y dependent variable
 * @param mu location
 * @param sigma scale
 * @param tau inverse scale
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if y is NaN, mu is infinite, sigma is negative or
 * infinite or tau is negative or infinite.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
          typename T_skewness_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_y_cl, T_loc_cl, T_scale_cl, T_skewness_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_scale_cl,
                                        T_skewness_cl>* = nullptr>
return_type_t<T_y_cl, T_loc_cl, T_scale_cl, T_skewness_cl>
skew_double_exponential_lpdf(const T_y_cl& y, const T_loc_cl& mu,
                             const T_scale_cl& sigma,
                             const T_skewness_cl& tau) {
  static const char* function = "skew_double_exponential_lpdf(OpenCL)";
  using T_partials_return
      = partials_return_t<T_y_cl, T_loc_cl, T_scale_cl, T_skewness_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma, "Inv_scale paramter",
                         tau);
  const size_t N = max_size(y, mu, sigma, tau);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_loc_cl, T_scale_cl,
                       T_skewness_cl>::value) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& mu_col = as_column_vector_or_scalar(mu);
  const auto& sigma_col = as_column_vector_or_scalar(sigma);
  const auto& tau_col = as_column_vector_or_scalar(tau);

  const auto& y_val = value_of(y_col);
  const auto& mu_val = value_of(mu_col);
  const auto& sigma_val = value_of(sigma_col);
  const auto& tau_val = value_of(tau_col);

  auto check_y_not_nan
      = check_cl(function, "Random variable", y_val, "not_nan");
  auto y_not_nan_expr = !isnan(y_val);
  auto check_mu_finite
      = check_cl(function, "Location parameter", mu_val, "finite");
  auto mu_finite_expr = isfinite(mu_val);
  auto check_sigma_positive_finite
      = check_cl(function, "Scale parameter", sigma_val, "positive finite");
  auto sigma_positive_finite_expr = isfinite(sigma_val) && sigma_val > 0;
  auto check_tau_bounded = check_cl(function, "Skewness parameter", tau_val,
                                    "in the interval [0, 1]");
  auto tau_bounded_expr = 0.0 < tau_val && tau_val < 1.0;

  auto inv_sigma = elt_divide(1.0, sigma_val);
  auto y_m_mu = y_val - mu_val;
  auto diff_sign = sign(y_m_mu);
  auto diff_sign_smaller_0 = diff_sign < 0.0;
  auto abs_diff_y_mu = fabs(y_m_mu);
  auto abs_diff_y_mu_over_sigma = elt_multiply(abs_diff_y_mu, inv_sigma);
  auto tmp = diff_sign_smaller_0 + elt_multiply(diff_sign, tau_val);
  auto expo = elt_multiply(tmp, abs_diff_y_mu_over_sigma);

  auto logp1 = -2.0 * expo;
  auto logp2 = static_select<include_summand<propto, T_scale_cl>::value>(
      logp1 - log(sigma_val), logp1);
  auto logp3 = static_select<include_summand<propto, T_skewness_cl>::value>(
      logp2 + log(tau_val) + log1m(tau_val), logp2);
  auto logp_expr = colwise_sum(logp3);

  auto mu_deriv = 2.0 * elt_multiply(tmp, elt_multiply(diff_sign, inv_sigma));
  auto y_deriv = -mu_deriv;
  auto sigma_deriv = -inv_sigma + 2.0 * elt_multiply(expo, inv_sigma);
  auto tau_deriv = elt_divide(1.0, tau_val) - elt_divide(1.0, 1.0 - tau_val)
                   - elt_multiply(diff_sign * 2.0, abs_diff_y_mu_over_sigma);

  matrix_cl<double> logp_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;
  matrix_cl<double> tau_deriv_cl;

  results(check_y_not_nan, check_mu_finite, check_sigma_positive_finite,
          check_tau_bounded, logp_cl, y_deriv_cl, mu_deriv_cl, sigma_deriv_cl,
          tau_deriv_cl)
      = expressions(y_not_nan_expr, mu_finite_expr, sigma_positive_finite_expr,
                    tau_bounded_expr, logp_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_loc_cl>::value>(mu_deriv),
                    calc_if<!is_constant<T_scale_cl>::value>(sigma_deriv),
                    calc_if<!is_constant<T_skewness_cl>::value>(tau_deriv));

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  if (include_summand<propto>::value) {
    logp += N * LOG_TWO;
  }

  auto ops_partials
      = make_partials_propagator(y_col, mu_col, sigma_col, tau_col);
  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  if (!is_constant<T_loc_cl>::value) {
    partials<1>(ops_partials) = std::move(mu_deriv_cl);
  }
  if (!is_constant<T_scale_cl>::value) {
    partials<2>(ops_partials) = std::move(sigma_deriv_cl);
  }
  if (!is_constant<T_skewness_cl>::value) {
    partials<3>(ops_partials) = std::move(tau_deriv_cl);
  }
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan

#endif
#endif
