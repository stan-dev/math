#ifndef STAN_MATH_OPENCL_PRIM_GUMBEL_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_GUMBEL_LPDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns the Gumbel log probability density for the given
 * location and scale. Given containers of matching sizes, returns the
 * log sum of densities.
 *
 * @tparam T_y_cl type of real parameter
 * @tparam T_loc_cl type of location parameter
 * @tparam T_scale_cl type of scale parameter
 * @param y real parameter
 * @param mu location parameter
 * @param beta scale parameter
 * @return log probability density or log sum of probability densities
 * @throw std::domain_error if y is nan, mu is infinite, or beta is nonpositive
 * @throw std::invalid_argument if container sizes mismatch
 */
template <
    bool propto, typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_loc_cl,
                                                T_scale_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_scale_cl>* = nullptr>
return_type_t<T_y_cl, T_loc_cl, T_scale_cl> gumbel_lpdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_scale_cl& beta) {
  using std::isfinite;
  using std::isnan;
  static const char* function = "gumbel_lpdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_loc_cl, T_scale_cl>;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", beta);
  const size_t N = max_size(y, mu, beta);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_loc_cl, T_scale_cl>::value) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& mu_col = as_column_vector_or_scalar(mu);
  const auto& beta_col = as_column_vector_or_scalar(beta);

  const auto& y_val = value_of(y_col);
  const auto& mu_val = value_of(mu_col);
  const auto& beta_val = value_of(beta_col);

  auto check_y_not_nan
      = check_cl(function, "Random variable", y_val, "not NaN");
  auto y_not_nan_expr = !isnan(y_val);
  auto check_mu_finite
      = check_cl(function, "Location parameter", mu_val, "finite");
  auto mu_finite_expr = isfinite(mu_val);
  auto check_beta_positive
      = check_cl(function, "Scale parameter", beta_val, "positive ");
  auto beta_positive_expr = beta_val > 0;

  auto inv_beta_expr = elt_divide(1.0, beta_val);
  auto y_minus_mu_over_beta_expr = elt_multiply(y_val - mu_val, inv_beta_expr);
  auto exp_y_m_mu_over_beta_expr = exp(-y_minus_mu_over_beta_expr);

  auto logp1_expr = -y_minus_mu_over_beta_expr - exp_y_m_mu_over_beta_expr;
  auto logp_expr
      = colwise_sum(static_select<include_summand<propto, T_scale_cl>::value>(
          logp1_expr - log(beta_val), logp1_expr));

  auto scaled_diff_expr
      = elt_multiply(inv_beta_expr, exp_y_m_mu_over_beta_expr) - inv_beta_expr;
  auto beta_deriv_expr
      = elt_multiply(-y_minus_mu_over_beta_expr, scaled_diff_expr)
        - inv_beta_expr;

  matrix_cl<double> logp_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> beta_deriv_cl;

  results(check_y_not_nan, check_mu_finite, check_beta_positive, logp_cl,
          y_deriv_cl, mu_deriv_cl, beta_deriv_cl)
      = expressions(y_not_nan_expr, mu_finite_expr, beta_positive_expr,
                    logp_expr,
                    calc_if<!is_constant<T_y_cl>::value>(scaled_diff_expr),
                    calc_if<!is_constant<T_loc_cl>::value>(-scaled_diff_expr),
                    calc_if<!is_constant<T_scale_cl>::value>(beta_deriv_expr));

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  auto ops_partials = make_partials_propagator(y_col, mu_col, beta_col);
  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  if (!is_constant<T_loc_cl>::value) {
    partials<1>(ops_partials) = std::move(mu_deriv_cl);
  }
  if (!is_constant<T_scale_cl>::value) {
    partials<2>(ops_partials) = std::move(beta_deriv_cl);
  }

  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
