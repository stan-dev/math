#ifndef STAN_MATH_OPENCL_PRIM_NORMAL_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_NORMAL_LPDF_HPP
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
 * The log of the normal density for the specified scalar(s) given
 * the specified mean(s) and deviation(s). y, mu, or sigma can
 * each be either a scalar or a vector matrix_cl. Any vector inputs
 * must be the same length.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each observation/mean/deviation triple.
 *
 * @tparam T_y_cl type of scalar
 * @tparam T_loc_cl type of location parameter
 * @tparam T_scale_cl type of scale parameter
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) location parameter(s)
 * for the normal distribution.
 * @param sigma (Sequence of) scale parameters for the normal distribution.
 * @return The log of the product of the densities.
 * @throw std::domain_error if the scale is not positive.
 */
template <
    bool propto, typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_loc_cl,
                                                T_scale_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_scale_cl>* = nullptr>
inline return_type_t<T_y_cl, T_loc_cl, T_scale_cl> normal_lpdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_scale_cl& sigma) {
  static const char* function = "normal_lpdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_loc_cl, T_scale_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
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

  auto check_y_not_nan
      = check_cl(function, "Random variable", y_val, "not NaN");
  auto y_not_nan = !isnan(y_val);
  auto check_mu_finite
      = check_cl(function, "Location parameter", mu_val, "finite");
  auto mu_finite = isfinite(mu_val);
  auto check_sigma_positive
      = check_cl(function, "Scale parameter", sigma_val, "positive");
  auto sigma_positive = 0 < sigma_val;

  auto inv_sigma = elt_divide(1., sigma_val);
  auto y_scaled = elt_multiply((y_val - mu_val), inv_sigma);
  auto y_scaled_sq = elt_multiply(y_scaled, y_scaled);

  auto logp1 = -0.5 * y_scaled_sq;
  auto logp_expr
      = colwise_sum(static_select<include_summand<propto, T_scale_cl>::value>(
          logp1 - log(sigma_val), logp1));

  auto scaled_diff = elt_multiply(inv_sigma, y_scaled);
  auto sigma_deriv = elt_multiply(inv_sigma, y_scaled_sq) - inv_sigma;

  matrix_cl<double> logp_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;

  results(check_y_not_nan, check_mu_finite, check_sigma_positive, logp_cl,
          y_deriv_cl, mu_deriv_cl, sigma_deriv_cl)
      = expressions(y_not_nan, mu_finite, sigma_positive, logp_expr,
                    calc_if<!is_constant<T_y_cl>::value>(-scaled_diff),
                    calc_if<!is_constant<T_loc_cl>::value>(scaled_diff),
                    calc_if<!is_constant<T_scale_cl>::value>(sigma_deriv));

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  if (include_summand<propto>::value) {
    logp += NEG_LOG_SQRT_TWO_PI * N;
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
