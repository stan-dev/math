#ifndef STAN_MATH_OPENCL_PRIM_WEIBULL_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_WEIBULL_LPDF_HPP
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
 * Returns the Weibull log probability density for the given
 * location and scale. Given containers of matching sizes, returns the
 * log sum of probability densities.
 *
 * @tparam T_y_cl type of real parameter
 * @tparam T_shape_cl type of shape parameter
 * @tparam T_scale_cl type of scale parameter
 * @param y real parameter
 * @param alpha shape parameter
 * @param sigma scale parameter
 * @return log probability density or log sum of probability densities
 * @throw std::domain_error if y is negative, alpha or sigma are nonpositive
 */
template <
    bool propto, typename T_y_cl, typename T_shape_cl, typename T_scale_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_shape_cl,
                                                T_scale_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_shape_cl, T_scale_cl>* = nullptr>
inline return_type_t<T_y_cl, T_shape_cl, T_scale_cl> weibull_lpdf(
    const T_y_cl& y, const T_shape_cl& alpha, const T_scale_cl& sigma) {
  static const char* function = "weibull_lpdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_shape_cl, T_scale_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Shape parameter",
                         alpha, "Scale parameter", sigma);
  const size_t N = max_size(y, alpha, sigma);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_shape_cl, T_scale_cl>::value) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& alpha_col = as_column_vector_or_scalar(alpha);
  const auto& sigma_col = as_column_vector_or_scalar(sigma);

  const auto& y_val = value_of(y_col);
  const auto& alpha_val = value_of(alpha_col);
  const auto& sigma_val = value_of(sigma_col);

  auto check_y_finite = check_cl(function, "Random variable", y_val, "finite");
  auto y_finite = isfinite(y_val);
  auto check_alpha_positive_finite
      = check_cl(function, "Shape parameter", alpha_val, "positive finite");
  auto alpha_positive_finite = isfinite(alpha_val) && 0 < alpha_val;
  auto check_sigma_positive_finite
      = check_cl(function, "Scale parameter", sigma_val, "positive finite");
  auto sigma_positive_finite = isfinite(sigma_val) && 0 < sigma_val;

  auto any_y_negative = colwise_max(cast<char>(y_val < 0));
  auto log_y = log(y_val);
  auto log_sigma = log(sigma_val);
  auto inv_sigma = elt_divide(1., sigma_val);
  auto y_div_sigma_pow_alpha = pow(elt_multiply(y_val, inv_sigma), alpha_val);

  auto logp1 = -y_div_sigma_pow_alpha;
  auto logp2 = static_select<include_summand<propto, T_shape_cl>::value>(
      logp1 + log(alpha_val), logp1);
  auto logp3
      = static_select<include_summand<propto, T_y_cl, T_shape_cl>::value>(
          logp2 + elt_multiply(alpha_val - 1.0, log_y), logp2);
  auto logp_expr = colwise_sum(
      static_select<include_summand<propto, T_shape_cl, T_scale_cl>::value>(
          logp3 - elt_multiply(alpha_val, log_sigma), logp3));

  auto one_m_y_div_sigma_pow_alpha = 1.0 - y_div_sigma_pow_alpha;
  auto y_deriv = elt_divide(
      elt_multiply(alpha_val, one_m_y_div_sigma_pow_alpha) - 1.0, y_val);
  auto alpha_deriv
      = elt_divide(1.0, alpha_val)
        + elt_multiply(one_m_y_div_sigma_pow_alpha, log_y - log_sigma);
  auto sigma_deriv = elt_multiply(elt_multiply(alpha_val, inv_sigma),
                                  -one_m_y_div_sigma_pow_alpha);

  matrix_cl<char> any_y_negative_cl;
  matrix_cl<double> logp_cl;
  matrix_cl<double> alpha_deriv_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;

  results(check_y_finite, check_alpha_positive_finite,
          check_sigma_positive_finite, any_y_negative_cl, logp_cl, y_deriv_cl,
          alpha_deriv_cl, sigma_deriv_cl)
      = expressions(y_finite, alpha_positive_finite, sigma_positive_finite,
                    any_y_negative, logp_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_shape_cl>::value>(alpha_deriv),
                    calc_if<!is_constant<T_scale_cl>::value>(sigma_deriv));

  if (from_matrix_cl(any_y_negative_cl).any()) {
    return LOG_ZERO;
  }

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  auto ops_partials = make_partials_propagator(y_col, alpha_col, sigma_col);

  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  if (!is_constant<T_shape_cl>::value) {
    partials<1>(ops_partials) = std::move(alpha_deriv_cl);
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
