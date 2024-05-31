#ifndef STAN_MATH_OPENCL_PRIM_DOUBLE_EXPONENTIAL_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_DOUBLE_EXPONENTIAL_LPDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sign.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/fun/sign.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns the double exponential log probability density function. Given
 * containers of matching sizes, returns the log sum of densities.
 *
 * @tparam T_y_cl type of real parameter.
 * @tparam T_loc_cl type of location parameter.
 * @tparam T_scale_cl type of scale parameter.
 * @param y real parameter
 * @param mu location parameter
 * @param sigma scale parameter
 * @return log probability density or log sum of probability densities
 * @throw std::domain_error if y is nan, mu is infinite, or sigma is nonpositive
 * @throw std::invalid_argument if container sizes mismatch
 */
template <
    bool propto, typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_loc_cl,
                                                T_scale_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_scale_cl>* = nullptr>
return_type_t<T_y_cl, T_loc_cl, T_scale_cl> double_exponential_lpdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_scale_cl& sigma) {
  static constexpr const char* function = "double_exponential_lpdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_loc_cl, T_scale_cl>;
  using std::isfinite;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Shape parameter", sigma);
  const size_t N = max_size(y, mu, sigma);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_loc_cl, T_scale_cl>::value) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& mu_col = as_column_vector_or_scalar(mu);
  const auto& sigma_col = as_column_vector_or_scalar(sigma);

  const auto& y_val = value_of(y_col);
  const auto& mu_val = value_of(mu_col);
  const auto& sigma_val = value_of(sigma_col);

  auto check_y_finite = check_cl(function, "Random variable", y_val, "finite");
  auto y_finite_expr = isfinite(y_val);
  auto check_mu_finite
      = check_cl(function, "Location parameter", mu_val, "finite");
  auto mu_finite_expr = isfinite(mu_val);
  auto check_sigma_positive_finite
      = check_cl(function, "Scale parameter", sigma_val, "positive finite");
  auto sigma_positive_finite_expr = sigma_val > 0 && isfinite(sigma_val);

  auto inv_sigma_expr = elt_divide(1.0, sigma_val);
  auto y_m_mu_expr = y_val - mu_val;
  auto abs_diff_y_mu_expr = fabs(y_m_mu_expr);
  auto scaled_diff_expr = elt_multiply(abs_diff_y_mu_expr, inv_sigma_expr);

  auto logp_expr
      = colwise_sum(-static_select<include_summand<propto, T_scale_cl>::value>(
          scaled_diff_expr + log(sigma_val), scaled_diff_expr));
  auto rep_deriv_expr = elt_multiply(sign(y_m_mu_expr), inv_sigma_expr);
  auto sigma_deriv_expr = elt_multiply(inv_sigma_expr, scaled_diff_expr - 1);

  matrix_cl<double> logp_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;
  results(check_y_finite, check_mu_finite, check_sigma_positive_finite, logp_cl,
          y_deriv_cl, mu_deriv_cl, sigma_deriv_cl)
      = expressions(y_finite_expr, mu_finite_expr, sigma_positive_finite_expr,
                    logp_expr, -rep_deriv_expr, rep_deriv_expr,
                    sigma_deriv_expr);

  T_partials_return logp = sum(from_matrix_cl(logp_cl));
  if (include_summand<propto>::value) {
    logp -= N * LOG_TWO;
  }
  auto ops_partials = make_partials_propagator(y_col, mu_col, sigma_col);

  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  if (!is_constant<T_loc_cl>::value) {
    partials<1>(ops_partials) = std::move(mu_deriv_cl);
  }
  if (!is_constant<T_scale_cl>::value) {
    partials<2>(ops_partials) = std::move(sigma_deriv_cl);
  }

  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan

#endif
#endif
