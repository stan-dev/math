#ifndef STAN_MATH_OPENCL_PRIM_DOUBLE_EXP_MOD_NORMAL_LCDF_HPP
#define STAN_MATH_OPENCL_PRIM_DOUBLE_EXP_MOD_NORMAL_LCDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/std_normal_lcdf.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns the exp mod normal log cumulative density
 * function. Given containers of matching sizes, returns the log sum of
 * probabilities.
 *
 * @tparam T_y_cl type of scalar outcome
 * @tparam T_loc_cl type of location
 * @tparam T_scale_cl type of scale
 * @tparam T_inv_scale_cl type of inverse scale
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) location(s).
 * @param sigma (Sequence of) scale(s).
 * @param lambda (Sequence of) inverse scale(s).
 * @return The log of the product of densities.
 */
template <typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
          typename T_inv_scale_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_y_cl, T_loc_cl, T_scale_cl, T_inv_scale_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_scale_cl,
                                        T_inv_scale_cl>* = nullptr>
inline return_type_t<T_y_cl, T_loc_cl, T_scale_cl, T_inv_scale_cl>
exp_mod_normal_lcdf(const T_y_cl& y, const T_loc_cl& mu,
                    const T_scale_cl& sigma, const T_inv_scale_cl& lambda) {
  static constexpr const char* function = "exp_mod_normal_lcdf(OpenCL)";
  using T_partials_return
      = partials_return_t<T_y_cl, T_loc_cl, T_scale_cl, T_inv_scale_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma, "Inv_scale parameter",
                         lambda);
  const size_t N = max_size(y, mu, sigma, lambda);
  if (N == 0) {
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
      = check_cl(function, "Random variable", y_val, "not NaN");
  auto y_not_nan_expr = !isnan(y_val);
  auto check_mu_finite
      = check_cl(function, "Location parameter", mu_val, "finite");
  auto mu_finite_expr = isfinite(mu_val);
  auto check_sigma_positive_finite
      = check_cl(function, "Scale parameter", sigma_val, "positive finite");
  auto sigma_positive_finite_expr = 0 < sigma_val && isfinite(sigma_val);
  auto check_lambda_positive_finite
      = check_cl(function, "Inv_cale parameter", lambda_val, "positive finite");
  auto lambda_positive_finite_expr = 0 < lambda_val && isfinite(lambda_val);

  auto any_y_neg_inf = colwise_max(cast<char>(y_val == NEGATIVE_INFTY));
  auto y_pos_inf = y_val == INFTY;
  auto sigma_inv = elt_divide(1.0, sigma_val);
  auto z = elt_multiply(y_val - mu_val, sigma_inv);
  auto a = elt_multiply(lambda_val, sigma_val);
  auto z_scaled = z * INV_SQRT_TWO;
  auto u_scaled = (z - a) * INV_SQRT_TWO;
  auto log_cdf_z = std_normal_lcdf_scaled_impl(z_scaled);
  auto log_cdf_u = std_normal_lcdf_scaled_impl(u_scaled);
  auto mills_z = std_normal_lcdf_dscaled_impl(z_scaled) * INV_SQRT_TWO;
  auto mills_u = std_normal_lcdf_dscaled_impl(u_scaled) * INV_SQRT_TWO;
  auto q = 0.5 * elt_multiply(a, a) - elt_multiply(a, z);
  auto erfc_arg = (a - z) * INV_SQRT_TWO;
  auto inv_two_erfc_arg_sq
      = elt_divide(0.5, elt_multiply(erfc_arg, erfc_arg));
  auto erfcx_series_5 = 105.0 - 945.0 * inv_two_erfc_arg_sq;
  auto erfcx_series_4
      = -15.0 + elt_multiply(inv_two_erfc_arg_sq, erfcx_series_5);
  auto erfcx_series_3
      = 3.0 + elt_multiply(inv_two_erfc_arg_sq, erfcx_series_4);
  auto erfcx_series_2
      = -1.0 + elt_multiply(inv_two_erfc_arg_sq, erfcx_series_3);
  auto erfcx_series
      = 1.0 + elt_multiply(inv_two_erfc_arg_sq, erfcx_series_2);
  auto erfcx_asymptotic
      = elt_divide(erfcx_series * INV_SQRT_PI, erfc_arg);
  auto erfcx_direct
      = elt_multiply(exp(elt_multiply(erfc_arg, erfc_arg)), erfc(erfc_arg));
  auto erfcx = select(erfc_arg >= 20.0, erfcx_asymptotic, erfcx_direct);
  auto use_erfcx = erfc_arg >= 5.0;
  auto stable_log_exp_cdf = select(
      use_erfcx,
      -0.5 * elt_multiply(z, z) + LOG_HALF + log(erfcx),
      q + log_cdf_u);
  auto inv_tail = elt_divide(1.0, a - z);
  auto inv_tail_sq = elt_multiply(inv_tail, inv_tail);
  auto mills_series_5 = 706.0 - 8162.0 * inv_tail_sq;
  auto mills_series_4
      = -74.0 + elt_multiply(inv_tail_sq, mills_series_5);
  auto mills_series_3
      = 10.0 + elt_multiply(inv_tail_sq, mills_series_4);
  auto mills_series_2
      = -2.0 + elt_multiply(inv_tail_sq, mills_series_3);
  auto mills_excess_asymptotic
      = elt_multiply(inv_tail,
                     1.0 + elt_multiply(inv_tail_sq, mills_series_2));
  auto mills_excess_erfcx
      = elt_divide(SQRT_TWO_OVER_SQRT_PI, erfcx) - (a - z);
  auto mills_excess
      = select(erfc_arg >= 20.0, mills_excess_asymptotic,
               select(use_erfcx, mills_excess_erfcx, mills_u - (a - z)));
  auto log_cdf_n = log_diff_exp(log_cdf_z, stable_log_exp_cdf);
  auto cdf_weight = exp(log_cdf_z - log_cdf_n);
  auto exp_cdf_weight = exp(stable_log_exp_cdf - log_cdf_n);
  auto dz_log_cdf
      = elt_multiply(cdf_weight, mills_z)
        + elt_multiply(exp_cdf_weight, z - mills_excess);
  auto da_log_cdf = elt_multiply(exp_cdf_weight, mills_excess);

  auto m0_factor = select(z < -4.0, mills_excess, z + mills_z);
  auto m1_over_m0 = elt_divide(
      0.5 * (elt_multiply(z, z) + 1.0 + elt_multiply(z, mills_z)),
      m0_factor);
  auto use_small_a = a < 1e-8 && m0_factor > 0.0
                     && fabs(elt_multiply(a, m1_over_m0)) < 1e-8;
  auto remainder = 1.0 - elt_multiply(a, m1_over_m0);
  auto small_log_cdf
      = log(a) + log_cdf_z + log(m0_factor) + log1p(-a * m1_over_m0);
  auto dm1_over_m0 = 1.0 - elt_divide(m1_over_m0, m0_factor);
  auto small_dz_log_cdf
      = elt_divide(1.0, m0_factor)
        - elt_divide(elt_multiply(a, dm1_over_m0), remainder);
  auto small_a_da_log_cdf
      = 1.0 - elt_divide(elt_multiply(a, m1_over_m0), remainder);
  auto small_da_log_cdf = elt_divide(small_a_da_log_cdf, a);
  auto stable_log_cdf = select(use_small_a, small_log_cdf, log_cdf_n);
  auto stable_dz_log_cdf
      = select(use_small_a, small_dz_log_cdf, dz_log_cdf);
  auto stable_da_log_cdf
      = select(use_small_a, small_da_log_cdf, da_log_cdf);
  auto stable_a_da_log_cdf
      = select(use_small_a, small_a_da_log_cdf, elt_multiply(a, da_log_cdf));

  auto cdf_log_expr = colwise_sum(select(y_pos_inf, 0.0, stable_log_cdf));
  auto y_deriv
      = select(y_pos_inf, 0.0, elt_multiply(stable_dz_log_cdf, sigma_inv));
  auto mu_deriv = -y_deriv;
  auto sigma_deriv = select(
      y_pos_inf, 0.0,
      elt_multiply(-elt_multiply(z, stable_dz_log_cdf)
                       + stable_a_da_log_cdf,
                   sigma_inv));
  auto lambda_deriv = select(
      y_pos_inf, 0.0, elt_multiply(sigma_val, stable_da_log_cdf));

  matrix_cl<char> any_y_neg_inf_cl;
  matrix_cl<double> cdf_log_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;
  matrix_cl<double> lambda_deriv_cl;

  results(check_y_not_nan, check_mu_finite, check_sigma_positive_finite,
          check_lambda_positive_finite, any_y_neg_inf_cl, cdf_log_cl,
          y_deriv_cl, mu_deriv_cl, sigma_deriv_cl, lambda_deriv_cl)
      = expressions(y_not_nan_expr, mu_finite_expr, sigma_positive_finite_expr,
                    lambda_positive_finite_expr, any_y_neg_inf, cdf_log_expr,
                    calc_if<is_autodiff_v<T_y_cl>>(y_deriv),
                    calc_if<is_autodiff_v<T_loc_cl>>(mu_deriv),
                    calc_if<is_autodiff_v<T_scale_cl>>(sigma_deriv),
                    calc_if<is_autodiff_v<T_inv_scale_cl>>(lambda_deriv));

  if (from_matrix_cl(any_y_neg_inf_cl).maxCoeff()) {
    return NEGATIVE_INFTY;
  }

  T_partials_return cdf_log = (from_matrix_cl(cdf_log_cl)).sum();

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
  return ops_partials.build(cdf_log);
}

}  // namespace math
}  // namespace stan
#endif
#endif
