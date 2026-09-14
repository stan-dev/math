#ifndef STAN_MATH_OPENCL_PRIM_EXP_MOD_NORMAL_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_EXP_MOD_NORMAL_LPDF_HPP
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
 * Returns the log PMF of the exp mod normal distribution. If
 * containers are supplied, returns the log sum of the probabilities.
 *
 * @tparam T_y_cl type of dependent variable
 * @tparam T_loc_cl type of location parameter
 * @tparam T_scale_cl type of scale parameter
 * @tparam T_inv_scale_cl type of inverse scale parameter
 * @param y dependent variable
 * @param mu location
 * @param sigma scale
 * @param lambda inverse scale
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if y is NaN, mu is infinite, sigma is negative or
 * infinite or lambda is negative or infinite.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
          typename T_inv_scale_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_y_cl, T_loc_cl, T_scale_cl, T_inv_scale_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_scale_cl,
                                        T_inv_scale_cl>* = nullptr>
inline return_type_t<T_y_cl, T_loc_cl, T_scale_cl, T_inv_scale_cl>
exp_mod_normal_lpdf(const T_y_cl& y, const T_loc_cl& mu,
                    const T_scale_cl& sigma, const T_inv_scale_cl& lambda) {
  static constexpr const char* function = "exp_mod_normal_lpdf(OpenCL)";
  using T_partials_return
      = partials_return_t<T_y_cl, T_loc_cl, T_scale_cl, T_inv_scale_cl>;
  using std::isfinite;
  using std::isinf;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma, "Inv_scale parameter",
                         lambda);
  const size_t N = max_size(y, mu, sigma, lambda);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_loc_cl, T_scale_cl,
                       T_inv_scale_cl>::value) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& mu_col = as_column_vector_or_scalar(mu);
  const auto& sigma_col = as_column_vector_or_scalar(sigma);
  const auto& lambda_col = as_column_vector_or_scalar(lambda);

  const auto& y_val = value_of(y_col);
  const auto& mu_val = value_of(mu_col);
  const auto& sigma_val = value_of(sigma_col);
  const auto& lambda_val = value_of(lambda_col);

  auto check_y_not_nan
      = check_cl(function, "Random variable", y_val, "not_nan");
  auto y_not_nan_expr = !isnan(y_val);
  auto check_mu_finite
      = check_cl(function, "Location parameter", mu_val, "finite");
  auto mu_finite_expr = isfinite(mu_val);
  auto check_sigma_positive_finite
      = check_cl(function, "Scale parameter", sigma_val, "positive finite");
  auto sigma_positive_finite_expr = isfinite(sigma_val) && sigma_val > 0;
  auto check_lambda_positive_finite = check_cl(function, "Inv_scale parameter",
                                               lambda_val, "positive finite");
  auto lambda_positive_finite_expr = isfinite(lambda_val) && lambda_val > 0;
  if constexpr (is_stan_scalar<T_y_cl>::value) {
    results(check_y_not_nan, check_mu_finite, check_sigma_positive_finite,
            check_lambda_positive_finite)
        = expressions(y_not_nan_expr, mu_finite_expr,
                      sigma_positive_finite_expr, lambda_positive_finite_expr);
    if (isinf(y_val)) {
      return NEGATIVE_INFTY;
    }
  } else {
    matrix_cl<char> any_y_inf_cl;
    results(check_y_not_nan, check_mu_finite, check_sigma_positive_finite,
            check_lambda_positive_finite, any_y_inf_cl)
        = expressions(y_not_nan_expr, mu_finite_expr,
                      sigma_positive_finite_expr, lambda_positive_finite_expr,
                      colwise_max(cast<char>(isinf(y_val))));
    if (from_matrix_cl(any_y_inf_cl).maxCoeff()) {
      return NEGATIVE_INFTY;
    }
  }

  auto inv_sigma_expr = elt_divide(1.0, sigma_val);
  auto z_expr = elt_multiply(y_val - mu_val, inv_sigma_expr);
  auto a_expr = elt_multiply(lambda_val, sigma_val);
  auto u_scaled_expr = (z_expr - a_expr) * INV_SQRT_TWO;
  auto log_cdf_u_expr = std_normal_lcdf_scaled_impl(u_scaled_expr);
  auto mills_u_expr
      = std_normal_lcdf_dscaled_impl(u_scaled_expr) * INV_SQRT_TWO;
  auto q_expr
      = 0.5 * elt_multiply(a_expr, a_expr) - elt_multiply(a_expr, z_expr);
  auto erfc_arg = (a_expr - z_expr) * INV_SQRT_TWO;
  auto inv_two_erfc_arg_sq = elt_divide(0.5, elt_multiply(erfc_arg, erfc_arg));
  auto erfcx_series
      = 1.0
        + elt_multiply(
            inv_two_erfc_arg_sq,
            -1.0
                + elt_multiply(
                    inv_two_erfc_arg_sq,
                    3.0
                        + elt_multiply(
                            inv_two_erfc_arg_sq,
                            -15.0
                                + elt_multiply(
                                    inv_two_erfc_arg_sq,
                                    105.0
                                        + elt_multiply(inv_two_erfc_arg_sq,
                                                       -945.0))))));
  auto erfcx_asymptotic = elt_divide(erfcx_series * INV_SQRT_PI, erfc_arg);
  auto erfcx_direct
      = elt_multiply(exp(elt_multiply(erfc_arg, erfc_arg)), erfc(erfc_arg));
  auto erfcx = select(erfc_arg >= 20.0, erfcx_asymptotic, erfcx_direct);
  auto use_erfcx = erfc_arg >= 5.0;
  auto log_exp_cdf_expr = select(
      use_erfcx, -0.5 * elt_multiply(z_expr, z_expr) + LOG_HALF + log(erfcx),
      q_expr + log_cdf_u_expr);
  auto inv_tail = elt_divide(1.0, a_expr - z_expr);
  auto inv_tail_sq = elt_multiply(inv_tail, inv_tail);
  auto mills_excess_asymptotic
      = elt_multiply(
          inv_tail,
          1.0
              + elt_multiply(
                  inv_tail_sq,
                  -2.0
                      + elt_multiply(
                          inv_tail_sq,
                          10.0
                              + elt_multiply(
                                  inv_tail_sq,
                                  -74.0
                                      + elt_multiply(inv_tail_sq,
                                                     706.0
                                                         - 8162.0
                                                               * inv_tail_sq))))));
  auto mills_excess_erfcx
      = elt_divide(SQRT_TWO_OVER_SQRT_PI, erfcx) - (a_expr - z_expr);
  auto mills_excess_expr = select(
      erfc_arg >= 20.0, mills_excess_asymptotic,
      select(use_erfcx, mills_excess_erfcx, mills_u_expr - (a_expr - z_expr)));
  auto logp1_expr = log_exp_cdf_expr + LOG_TWO;
  auto logp_expr = colwise_sum(
      static_select<include_summand<propto, T_inv_scale_cl>::value>(
          logp1_expr + log(lambda_val), logp1_expr));

  auto dz_expr = -z_expr + mills_excess_expr;
  auto da_expr = -mills_excess_expr;
  auto deriv_sigma_expr = elt_multiply(
      -elt_multiply(z_expr, dz_expr) + elt_multiply(a_expr, da_expr),
      inv_sigma_expr);
  auto deriv_lambda_expr
      = elt_divide(1.0, lambda_val) + elt_multiply(sigma_val, da_expr);

  matrix_cl<double> logp_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;
  matrix_cl<double> lambda_deriv_cl;

  results(logp_cl, y_deriv_cl, mu_deriv_cl, sigma_deriv_cl, lambda_deriv_cl)
      = expressions(
          logp_expr,
          calc_if<is_autodiff_v<T_y_cl>>(elt_multiply(dz_expr, inv_sigma_expr)),
          calc_if<is_autodiff_v<T_loc_cl>>(
              -elt_multiply(dz_expr, inv_sigma_expr)),
          calc_if<is_autodiff_v<T_scale_cl>>(deriv_sigma_expr),
          calc_if<is_autodiff_v<T_inv_scale_cl>>(deriv_lambda_expr));

  T_partials_return logp = sum(from_matrix_cl(logp_cl));
  if constexpr (include_summand<propto>::value) {
    logp -= LOG_TWO * N;
  }

  auto ops_partials
      = make_partials_propagator(y_col, mu_col, sigma_col, lambda_col);
  if constexpr (is_autodiff_v<T_y_cl>) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  if constexpr (is_autodiff_v<T_loc_cl>) {
    partials<1>(ops_partials) = std::move(mu_deriv_cl);
  }
  if constexpr (is_autodiff_v<T_scale_cl>) {
    partials<2>(ops_partials) = std::move(sigma_deriv_cl);
  }
  if constexpr (is_autodiff_v<T_inv_scale_cl>) {
    partials<3>(ops_partials) = std::move(lambda_deriv_cl);
  }
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan

#endif
#endif
