#ifndef STAN_MATH_OPENCL_PRIM_SKEW_NORMAL_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_SKEW_NORMAL_LPDF_HPP
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
 * The log of the skew normal density for the specified scalar(s) given
 * the specified mean(s), deviation(s) and shape(s). y, mu, sigma, or alpha can
 * each be either a scalar or a vector matrix_cl. Any vector inputs
 * must be the same length.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each observation/mean/deviation quadruple.
 *
 * @tparam T_y_cl type of scalar
 * @tparam T_loc_cl type of location parameter
 * @tparam T_scale_cl type of scale parameter
 * @tparam T_shape_cl type of shape parameter
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) location parameter(s)
 * @param sigma (Sequence of) scale parameter(s)
 * @param alpha (Sequence of) shape parameter(s)
 * @return The log of the product of the densities.
 * @throw std::domain_error if the scale is not positive.
 */
template <bool propto, typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
          typename T_shape_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_y_cl, T_loc_cl, T_scale_cl, T_shape_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_scale_cl,
                                        T_shape_cl>* = nullptr>
inline return_type_t<T_y_cl, T_loc_cl, T_scale_cl, T_shape_cl> skew_normal_lpdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_scale_cl& sigma,
    const T_shape_cl& alpha) {
  static const char* function = "skew_normal_lpdf(OpenCL)";
  using T_partials_return
      = partials_return_t<T_y_cl, T_loc_cl, T_scale_cl, T_shape_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma, "Shape paramter", alpha);
  const size_t N = max_size(y, mu, sigma, alpha);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_loc_cl, T_scale_cl,
                       T_shape_cl>::value) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& mu_col = as_column_vector_or_scalar(mu);
  const auto& sigma_col = as_column_vector_or_scalar(sigma);
  const auto& alpha_col = as_column_vector_or_scalar(alpha);

  const auto& y_val = value_of(y_col);
  const auto& mu_val = value_of(mu_col);
  const auto& sigma_val = value_of(sigma_col);
  const auto& alpha_val = value_of(alpha_col);

  auto check_y_not_nan
      = check_cl(function, "Random variable", y_val, "not NaN");
  auto y_not_nan = !isnan(y_val);
  auto check_mu_finite
      = check_cl(function, "Location parameter", mu_val, "finite");
  auto mu_finite = isfinite(mu_val);
  auto check_sigma_positive
      = check_cl(function, "Scale parameter", sigma_val, "positive");
  auto sigma_positive = 0 < sigma_val;
  auto check_alpha_finite
      = check_cl(function, "Shape parameter", alpha_val, "finite");
  auto alpha_finite = isfinite(alpha_val);

  auto inv_sigma = elt_divide(1., sigma_val);
  auto y_minus_mu_over_sigma = elt_multiply((y_val - mu_val), inv_sigma);
  auto log_erfc_alpha_z = log(
      erfc(elt_multiply(alpha_val, y_minus_mu_over_sigma) * -INV_SQRT_TWO));

  auto logp1 = log_erfc_alpha_z;
  auto logp2 = static_select<include_summand<propto, T_scale_cl>::value>(
      logp1 - log(sigma_val), logp1);
  auto logp_expr = colwise_sum(
      static_select<
          include_summand<propto, T_y_cl, T_loc_cl, T_scale_cl>::value>(
          logp2
              - elt_multiply(y_minus_mu_over_sigma, y_minus_mu_over_sigma)
                    * 0.5,
          logp2));

  auto scaled = elt_multiply(alpha_val, y_minus_mu_over_sigma) * INV_SQRT_TWO;
  auto deriv_logerf = SQRT_TWO_OVER_SQRT_PI
                      * exp(-elt_multiply(scaled, scaled) - log_erfc_alpha_z);
  auto y_loc_deriv = elt_multiply(
      y_minus_mu_over_sigma - elt_multiply(deriv_logerf, alpha_val), inv_sigma);
  auto sigma_deriv
      = elt_multiply(elt_multiply(y_minus_mu_over_sigma
                                      - elt_multiply(deriv_logerf, alpha_val),
                                  y_minus_mu_over_sigma)
                         - 1,
                     inv_sigma);
  auto alpha_deriv = elt_multiply(deriv_logerf, y_minus_mu_over_sigma);

  matrix_cl<double> logp_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;
  matrix_cl<double> alpha_deriv_cl;

  results(check_y_not_nan, check_mu_finite, check_sigma_positive,
          check_alpha_finite, logp_cl, y_deriv_cl, mu_deriv_cl, sigma_deriv_cl,
          alpha_deriv_cl)
      = expressions(y_not_nan, mu_finite, sigma_positive, alpha_finite,
                    logp_expr,
                    calc_if<!is_constant<T_y_cl>::value>(-y_loc_deriv),
                    calc_if<!is_constant<T_loc_cl>::value>(y_loc_deriv),
                    calc_if<!is_constant<T_scale_cl>::value>(sigma_deriv),
                    calc_if<!is_constant<T_shape_cl>::value>(alpha_deriv));

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  if (include_summand<propto>::value) {
    logp -= HALF_LOG_TWO_PI * N;
  }

  operands_and_partials<decltype(y_col), decltype(mu_col), decltype(sigma_col),
                        decltype(alpha_col)>
      ops_partials(y_col, mu_col, sigma_col, alpha_col);

  if (!is_constant<T_y_cl>::value) {
    ops_partials.edge1_.partials_ = std::move(y_deriv_cl);
  }
  if (!is_constant<T_loc_cl>::value) {
    ops_partials.edge2_.partials_ = std::move(mu_deriv_cl);
  }
  if (!is_constant<T_scale_cl>::value) {
    ops_partials.edge3_.partials_ = std::move(sigma_deriv_cl);
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
