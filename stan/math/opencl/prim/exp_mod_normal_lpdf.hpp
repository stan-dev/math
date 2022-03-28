#ifndef STAN_MATH_OPENCL_PRIM_EXP_MOD_NORMAL_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_EXP_MOD_NORMAL_LPDF_HPP
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
return_type_t<T_y_cl, T_loc_cl, T_scale_cl, T_inv_scale_cl> exp_mod_normal_lpdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_scale_cl& sigma,
    const T_inv_scale_cl& lambda) {
  static const char* function = "exp_mod_normal_lpdf(OpenCL)";
  using T_partials_return
      = partials_return_t<T_y_cl, T_loc_cl, T_scale_cl, T_inv_scale_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma, "Inv_scale paramter",
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

  auto inv_sigma_expr = elt_divide(1.0, sigma_val);
  auto sigma_sq_expr = elt_multiply(sigma_val, sigma_val);
  auto lambda_sigma_sq_expr = elt_multiply(lambda_val, sigma_sq_expr);
  auto mu_minus_y_expr = mu_val - y_val;
  auto inner_term_expr = elt_multiply(mu_minus_y_expr + lambda_sigma_sq_expr,
                                      INV_SQRT_TWO * inv_sigma_expr);
  auto erfc_calc_expr = erfc(inner_term_expr);
  auto logp1_expr
      = elt_multiply(lambda_val, mu_minus_y_expr + 0.5 * lambda_sigma_sq_expr)
        + log(erfc_calc_expr);
  auto logp_expr = colwise_sum(
      static_select<include_summand<propto, T_inv_scale_cl>::value>(
          logp1_expr + log(lambda_val), logp1_expr));

  auto deriv_logerfc_expr
      = elt_divide(-SQRT_TWO_OVER_SQRT_PI
                       * exp(-elt_multiply(inner_term_expr, inner_term_expr)),
                   erfc_calc_expr);
  auto deriv_expr
      = lambda_val + elt_multiply(deriv_logerfc_expr, inv_sigma_expr);
  auto deriv_sigma_expr
      = elt_multiply(sigma_val, elt_multiply(lambda_val, lambda_val))
        + elt_multiply(
            deriv_logerfc_expr,
            (lambda_val - elt_divide(mu_minus_y_expr, sigma_sq_expr)));
  auto deriv_lambda_expr = elt_divide(1.0, lambda_val) + lambda_sigma_sq_expr
                           + mu_minus_y_expr
                           + elt_multiply(deriv_logerfc_expr, sigma_val);

  matrix_cl<double> logp_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;
  matrix_cl<double> lambda_deriv_cl;

  results(check_y_not_nan, check_mu_finite, check_sigma_positive_finite,
          check_lambda_positive_finite, logp_cl, y_deriv_cl, mu_deriv_cl,
          sigma_deriv_cl, lambda_deriv_cl)
      = expressions(
          y_not_nan_expr, mu_finite_expr, sigma_positive_finite_expr,
          lambda_positive_finite_expr, logp_expr,
          calc_if<!is_constant<T_y_cl>::value>(-deriv_expr),
          calc_if<!is_constant<T_loc_cl>::value>(deriv_expr),
          calc_if<!is_constant<T_scale_cl>::value>(deriv_sigma_expr),
          calc_if<!is_constant<T_inv_scale_cl>::value>(deriv_lambda_expr));

  T_partials_return logp = sum(from_matrix_cl(logp_cl));
  if (include_summand<propto>::value) {
    logp -= LOG_TWO * N;
  }

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
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan

#endif
#endif
