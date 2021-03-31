#ifndef STAN_MATH_OPENCL_PRIM_DOUBLE_EXP_MOD_NORMAL_LCCDF_HPP
#define STAN_MATH_OPENCL_PRIM_DOUBLE_EXP_MOD_NORMAL_LCCDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns the exp mod normal log complementary cumulative density
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
return_type_t<T_y_cl, T_loc_cl, T_scale_cl, T_inv_scale_cl>
exp_mod_normal_lccdf(const T_y_cl& y, const T_loc_cl& mu,
                     const T_scale_cl& sigma, const T_inv_scale_cl& lambda) {
  static const char* function = "exp_mod_normal_lccdf(OpenCL)";
  using T_partials_return
      = partials_return_t<T_y_cl, T_loc_cl, T_scale_cl, T_inv_scale_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
  const size_t N = max_size(y, mu, sigma);
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

  auto any_y_neg_inf
      = colwise_max(constant(0, N, 1) + (y_val == NEGATIVE_INFTY));
  auto any_y_pos_inf = colwise_max(constant(0, N, 1) + (y_val == INFTY));
  auto inv_sigma = elt_divide(1.0, sigma_val);
  auto diff = y_val - mu_val;
  auto scaled_diff = elt_multiply(diff, inv_sigma * INV_SQRT_TWO);
  auto v = elt_multiply(lambda_val, sigma_val);
  auto scaled_diff_diff = scaled_diff - v * INV_SQRT_TWO;
  auto erf_calc = 0.5 * (1.0 + erf(scaled_diff_diff));
  auto exp_term = exp(0.5 * square(v) - elt_multiply(lambda_val, diff));
  auto ccdf_n = 0.5 - 0.5 * erf(scaled_diff) + elt_multiply(exp_term, erf_calc);
  auto ccdf_log_expr = colwise_sum(log(ccdf_n));

  auto exp_term_2 = exp(-square(scaled_diff_diff));
  auto deriv_1 = elt_multiply(elt_multiply(lambda_val, exp_term), erf_calc);
  auto deriv_2 = INV_SQRT_TWO_PI
                 * elt_multiply(elt_multiply(exp_term, exp_term_2), inv_sigma);
  auto deriv_3
      = INV_SQRT_TWO_PI * elt_multiply(exp(-square(scaled_diff)), inv_sigma);
  auto mu_deriv = elt_divide(deriv_1 - deriv_2 + deriv_3, ccdf_n);
  auto y_deriv = -mu_deriv;
  auto sigma_deriv = elt_divide(
      elt_multiply(deriv_1 - deriv_2, v)
          + elt_multiply(deriv_3 - deriv_2, scaled_diff) * SQRT_TWO,
      ccdf_n);
  auto lambda_deriv = elt_divide(
      elt_multiply(exp_term,
                   elt_multiply(elt_multiply(v, sigma_val) - diff, erf_calc)
                       - INV_SQRT_TWO_PI * elt_multiply(sigma_val, exp_term_2)),
      ccdf_n);

  matrix_cl<double> any_y_neg_inf_cl;
  matrix_cl<double> any_y_pos_inf_cl;
  matrix_cl<double> ccdf_log_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;
  matrix_cl<double> lambda_deriv_cl;

  results(check_y_not_nan, check_mu_finite, check_sigma_positive_finite,
          check_lambda_positive_finite, any_y_neg_inf_cl, any_y_pos_inf_cl,
          ccdf_log_cl, y_deriv_cl, mu_deriv_cl, sigma_deriv_cl, lambda_deriv_cl)
      = expressions(y_not_nan_expr, mu_finite_expr, sigma_positive_finite_expr,
                    lambda_positive_finite_expr, any_y_neg_inf, any_y_pos_inf,
                    ccdf_log_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_loc_cl>::value>(mu_deriv),
                    calc_if<!is_constant<T_scale_cl>::value>(sigma_deriv),
                    calc_if<!is_constant<T_inv_scale_cl>::value>(lambda_deriv));

  if (from_matrix_cl(any_y_pos_inf_cl).maxCoeff()) {
    return NEGATIVE_INFTY;
  }

  if (from_matrix_cl(any_y_neg_inf_cl).maxCoeff()) {
    return 0.0;
  }

  T_partials_return ccdf_log = (from_matrix_cl(ccdf_log_cl)).sum();

  operands_and_partials<decltype(y_col), decltype(mu_col), decltype(sigma_col),
                        decltype(lambda_col)>
      ops_partials(y_col, mu_col, sigma_col, lambda_col);

  if (!is_constant<T_y_cl>::value) {
    ops_partials.edge1_.partials_ = std::move(y_deriv_cl);
  }
  if (!is_constant<T_loc_cl>::value) {
    ops_partials.edge2_.partials_ = std::move(mu_deriv_cl);
  }
  if (!is_constant<T_scale_cl>::value) {
    ops_partials.edge3_.partials_ = std::move(sigma_deriv_cl);
  }
  if (!is_constant<T_inv_scale_cl>::value) {
    ops_partials.edge4_.partials_ = std::move(lambda_deriv_cl);
  }
  return ops_partials.build(ccdf_log);
}

}  // namespace math
}  // namespace stan
#endif
#endif