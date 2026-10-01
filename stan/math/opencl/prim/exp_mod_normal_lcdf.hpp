#ifndef STAN_MATH_OPENCL_PRIM_DOUBLE_EXP_MOD_NORMAL_LCDF_HPP
#define STAN_MATH_OPENCL_PRIM_DOUBLE_EXP_MOD_NORMAL_LCDF_HPP
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
namespace internal {

/** Same log-space formulation as the prim exp_mod_normal_lcdf_impl. */
template <bool upper, typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
          typename T_inv_scale_cl>
inline return_type_t<T_y_cl, T_loc_cl, T_scale_cl, T_inv_scale_cl>
exp_mod_normal_lcdf_opencl_impl(const char* function, const T_y_cl& y,
                                const T_loc_cl& mu, const T_scale_cl& sigma,
                                const T_inv_scale_cl& lambda) {
  constexpr double sign = upper ? -1.0 : 1.0;
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
  auto check_lambda_positive_finite = check_cl(function, "Inv_scale parameter",
                                               lambda_val, "positive finite");
  auto lambda_positive_finite_expr = 0 < lambda_val && isfinite(lambda_val);

  auto any_y_neg_inf = colwise_max(cast<char>(y_val == NEGATIVE_INFTY));
  auto any_y_pos_inf = colwise_max(cast<char>(y_val == INFTY));
  auto diff = y_val - mu_val;
  auto z = elt_divide(diff, sigma_val);
  auto v = elt_multiply(lambda_val, sigma_val);
  auto z_a = sign * z;
  auto z_b = z - v;
  auto log_a = math::std_normal_lcdf_impl(z_a);
  auto log_b = 0.5 * square(v) - elt_multiply(lambda_val, diff)
               + math::std_normal_lcdf_impl(z_b);
  auto lp = [&]() {
    if constexpr (upper) {
      return fmax(log_a, log_b) + log1p_exp(-fabs(log_a - log_b));
    } else {
      return log_diff_exp(log_a, log_b);
    }
  }();
  auto cdf_log_expr = colwise_sum(lp);

  auto w_b = -sign * exp(log_b - lp);
  auto s_b = elt_multiply(w_b, std_normal_lcdf_derivative(z_b));
  auto s = elt_divide(
      sign * elt_multiply(exp(log_a - lp), std_normal_lcdf_derivative(z_a))
          + s_b,
      sigma_val);
  auto q = elt_multiply(w_b, v) - s_b;
  auto y_deriv = s - elt_multiply(w_b, lambda_val);
  auto mu_deriv = -y_deriv;
  auto sigma_deriv
      = select(s == 0, 0.0, -elt_multiply(s, z)) + elt_multiply(lambda_val, q);
  auto lambda_deriv = elt_multiply(sigma_val, q) - elt_multiply(w_b, diff);

  matrix_cl<char> any_y_neg_inf_cl;
  matrix_cl<char> any_y_pos_inf_cl;
  matrix_cl<double> cdf_log_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;
  matrix_cl<double> lambda_deriv_cl;

  results(check_y_not_nan, check_mu_finite, check_sigma_positive_finite,
          check_lambda_positive_finite, any_y_neg_inf_cl, any_y_pos_inf_cl,
          cdf_log_cl, y_deriv_cl, mu_deriv_cl, sigma_deriv_cl, lambda_deriv_cl)
      = expressions(y_not_nan_expr, mu_finite_expr, sigma_positive_finite_expr,
                    lambda_positive_finite_expr, any_y_neg_inf, any_y_pos_inf,
                    cdf_log_expr, calc_if<is_autodiff_v<T_y_cl>>(y_deriv),
                    calc_if<is_autodiff_v<T_loc_cl>>(mu_deriv),
                    calc_if<is_autodiff_v<T_scale_cl>>(sigma_deriv),
                    calc_if<is_autodiff_v<T_inv_scale_cl>>(lambda_deriv));

  if (from_matrix_cl(any_y_neg_inf_cl).maxCoeff()) {
    return upper ? 0.0 : NEGATIVE_INFTY;
  }
  if (from_matrix_cl(any_y_pos_inf_cl).maxCoeff()) {
    return upper ? NEGATIVE_INFTY : 0.0;
  }

  double cdf_log = sum(from_matrix_cl(cdf_log_cl));

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

}  // namespace internal

template <typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
          typename T_inv_scale_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_y_cl, T_loc_cl, T_scale_cl, T_inv_scale_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_scale_cl,
                                        T_inv_scale_cl>* = nullptr>
inline return_type_t<T_y_cl, T_loc_cl, T_scale_cl, T_inv_scale_cl>
exp_mod_normal_lcdf(const T_y_cl& y, const T_loc_cl& mu,
                    const T_scale_cl& sigma, const T_inv_scale_cl& lambda) {
  return internal::exp_mod_normal_lcdf_opencl_impl<false>(
      "exp_mod_normal_lcdf(OpenCL)", y, mu, sigma, lambda);
}

}  // namespace math
}  // namespace stan
#endif
#endif
