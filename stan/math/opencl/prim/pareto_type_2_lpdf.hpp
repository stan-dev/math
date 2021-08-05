#ifndef STAN_MATH_OPENCL_PRIM_PARETO_TYPE_2_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_PARETO_TYPE_2_LPDF_HPP
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
 * Returns the log PMF of the Pareto type 2 distribution. If
 * containers are supplied, returns the log sum of the probabilities.
 *
 * @tparam T_y_cl type of dependent variable
 * @tparam T_loc_cl type of location parameter
 * @tparam T_scale_cl type of scale parameter
 * @tparam T_shape_cl type of inverse scale parameter
 * @param y dependent variable
 * @param mu location
 * @param lambda scale
 * @param alpha inverse scale
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if y is NaN, mu is infinite, lambda is negative or
 * infinite or alpha is negative or infinite.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
          typename T_shape_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_y_cl, T_loc_cl, T_scale_cl, T_shape_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_scale_cl,
                                        T_shape_cl>* = nullptr>
return_type_t<T_y_cl, T_loc_cl, T_scale_cl, T_shape_cl> pareto_type_2_lpdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_scale_cl& lambda,
    const T_shape_cl& alpha) {
  static const char* function = "pareto_type_2_lpdf(OpenCL)";
  using T_partials_return
      = partials_return_t<T_y_cl, T_loc_cl, T_scale_cl, T_shape_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", alpha, "Shape parameter",
                         alpha);
  const size_t N = max_size(y, mu, lambda, alpha);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_loc_cl, T_scale_cl,
                       T_shape_cl>::value) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& mu_col = as_column_vector_or_scalar(mu);
  const auto& lambda_col = as_column_vector_or_scalar(lambda);
  const auto& alpha_col = as_column_vector_or_scalar(alpha);

  const auto& y_val = value_of(y_col);
  const auto& mu_val = value_of(mu_col);
  const auto& lambda_val = value_of(lambda_col);
  const auto& alpha_val = value_of(alpha_col);

  auto y_minus_mu = y_val - mu_val;
  auto check_y_ge_mu
      = check_cl(function, "Random variable minus location parameter",
                 y_minus_mu, "greater or equal than zero");
  auto y_ge_mu = y_minus_mu >= 0;
  auto check_lambda_positive_finite
      = check_cl(function, "Scale parameter", lambda_val, "positive finite");
  auto lambda_positive_finite = isfinite(lambda_val) && lambda_val > 0;
  auto check_alpha_positive_finite
      = check_cl(function, "Shape parameter", alpha_val, "positive finite");
  auto alpha_positive_finite = isfinite(alpha_val) && alpha_val > 0;

  auto log1p_scaled_diff = log1p(elt_divide(y_minus_mu, lambda_val));

  auto logp1 = static_select<include_summand<propto, T_shape_cl>::value>(
      log(alpha_val), constant(0, N, 1));
  auto logp2 = static_select<include_summand<propto, T_scale_cl>::value>(
      logp1 - log(lambda_val), logp1);
  auto logp_expr = colwise_sum(
      static_select<include_summand<propto, T_y_cl, T_loc_cl, T_scale_cl,
                                    T_shape_cl>::value>(
          logp2 - elt_multiply(alpha_val + 1.0, log1p_scaled_diff), logp2));

  auto inv_sum = elt_divide(1.0, lambda_val + y_minus_mu);
  auto alpha_div_sum = elt_multiply(alpha_val, inv_sum);

  auto deriv_y_mu = inv_sum + alpha_div_sum;
  auto deriv_lambda
      = elt_divide(elt_multiply(alpha_div_sum, y_minus_mu), lambda_val)
        - inv_sum;
  auto deriv_alpha = elt_divide(1.0, alpha_val) - log1p_scaled_diff;

  matrix_cl<double> logp_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> lambda_deriv_cl;
  matrix_cl<double> alpha_deriv_cl;

  results(check_y_ge_mu, check_lambda_positive_finite,
          check_alpha_positive_finite, logp_cl, y_deriv_cl, mu_deriv_cl,
          lambda_deriv_cl, alpha_deriv_cl)
      = expressions(y_ge_mu, lambda_positive_finite, alpha_positive_finite,
                    logp_expr,
                    calc_if<!is_constant<T_y_cl>::value>(-deriv_y_mu),
                    calc_if<!is_constant<T_loc_cl>::value>(deriv_y_mu),
                    calc_if<!is_constant<T_scale_cl>::value>(deriv_lambda),
                    calc_if<!is_constant<T_shape_cl>::value>(deriv_alpha));

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  operands_and_partials<decltype(y_col), decltype(mu_col), decltype(lambda_col),
                        decltype(alpha_col)>
      ops_partials(y_col, mu_col, lambda_col, alpha_col);
  if (!is_constant<T_y_cl>::value) {
    ops_partials.edge1_.partials_ = std::move(y_deriv_cl);
  }
  if (!is_constant<T_loc_cl>::value) {
    ops_partials.edge2_.partials_ = std::move(mu_deriv_cl);
  }
  if (!is_constant<T_scale_cl>::value) {
    ops_partials.edge3_.partials_ = std::move(lambda_deriv_cl);
  }
  if (!is_constant<T_shape_cl>::value) {
    ops_partials.edge4_.partials_ = std::move(alpha_deriv_cl);
  }
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan

#endif
#endif
