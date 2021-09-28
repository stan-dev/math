#ifndef STAN_MATH_OPENCL_PRIM_DOUBLE_PARETO_TYPE_2_CDF_HPP
#define STAN_MATH_OPENCL_PRIM_DOUBLE_PARETO_TYPE_2_CDF_HPP
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
 * Returns the pareto type 2 cumulative density function. Given
 * containers of matching sizes, returns the product of probabilities.
 *
 * @tparam T_y_cl type of scalar outcome
 * @tparam T_loc_cl type of location
 * @tparam T_scale_cl type of scale
 * @tparam T_shape_cl type of inverse scale
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) location(s).
 * @param lambda (Sequence of) scale(s).
 * @param alpha (Sequence of) shape(s).
 * @return The log of the product of densities.
 */
template <typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
          typename T_shape_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_y_cl, T_loc_cl, T_scale_cl, T_shape_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_scale_cl,
                                        T_shape_cl>* = nullptr>
return_type_t<T_y_cl, T_loc_cl, T_scale_cl, T_shape_cl> pareto_type_2_cdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_scale_cl& lambda,
    const T_shape_cl& alpha) {
  static const char* function = "pareto_type_2_cdf(OpenCL)";
  using T_partials_return
      = partials_return_t<T_y_cl, T_loc_cl, T_scale_cl, T_shape_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", lambda, "Shape parameter",
                         alpha);
  const size_t N = max_size(y, mu, lambda, alpha);
  if (N == 0) {
    return 1.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& mu_col = as_column_vector_or_scalar(mu);
  const auto& lambda_col = as_column_vector_or_scalar(lambda);
  const auto& alpha_col = as_column_vector_or_scalar(alpha);

  const auto& y_val = value_of(y_col);
  const auto& mu_val = value_of(mu_col);
  const auto& lambda_val = value_of(lambda_col);
  const auto& alpha_val = value_of(alpha_col);

  auto check_y_nonnegative
      = check_cl(function, "Random variable", y_val, "nonnegative");
  auto y_nonnegative_expr = 0 <= y_val;
  auto check_lambda_positive_finite
      = check_cl(function, "Scale parameter", lambda_val, "positive finite");
  auto lambda_positive_finite_expr = 0 < lambda_val && isfinite(lambda_val);
  auto check_alpha_positive_finite
      = check_cl(function, "Shape parameter", alpha_val, "positive finite");
  auto alpha_positive_finite_expr = 0 < alpha_val && isfinite(alpha_val);
  auto diff = y_val - mu_val;
  auto check_diff_nonnegative
      = check_cl(function, "Random variable minus location parameter", diff,
                 "nonnegative");
  auto diff_nonnegative_expr = 0 <= diff;

  auto summed = lambda_val + diff;
  auto temp = elt_divide(summed, lambda_val);
  auto p1_pow_alpha = pow(temp, -alpha_val);
  auto cdf_expr = colwise_prod(1.0 - p1_pow_alpha);

  auto inv_cdf_n = elt_divide(1.0, 1.0 - p1_pow_alpha);
  auto y_deriv1 = elt_multiply(elt_divide(p1_pow_alpha, summed),
                               elt_multiply(alpha_val, inv_cdf_n));
  auto lambda_deriv1 = elt_multiply(elt_divide(diff, -lambda_val), y_deriv1);
  auto alpha_deriv1
      = elt_multiply(elt_multiply(log(temp), p1_pow_alpha), inv_cdf_n);

  matrix_cl<double> cdf_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> lambda_deriv_cl;
  matrix_cl<double> alpha_deriv_cl;

  results(check_y_nonnegative, check_lambda_positive_finite,
          check_alpha_positive_finite, check_diff_nonnegative, cdf_cl,
          mu_deriv_cl, lambda_deriv_cl, alpha_deriv_cl)
      = expressions(
          y_nonnegative_expr, lambda_positive_finite_expr,
          alpha_positive_finite_expr, diff_nonnegative_expr, cdf_expr,
          calc_if<!is_constant_all<T_y_cl, T_loc_cl>::value>(y_deriv1),
          calc_if<!is_constant<T_scale_cl>::value>(lambda_deriv1),
          calc_if<!is_constant<T_shape_cl>::value>(alpha_deriv1));

  T_partials_return cdf = (from_matrix_cl(cdf_cl)).prod();

  operands_and_partials<decltype(y_col), decltype(mu_col), decltype(lambda_col),
                        decltype(alpha_col)>
      ops_partials(y_col, mu_col, lambda_col, alpha_col);
  if (!is_constant_all<T_y_cl, T_loc_cl, T_scale_cl, T_shape_cl>::value) {
    auto y_deriv = mu_deriv_cl * cdf;
    auto mu_deriv = -y_deriv;
    auto lambda_deriv = lambda_deriv_cl * cdf;
    auto alpha_deriv = alpha_deriv_cl * cdf;

    results(y_deriv_cl, mu_deriv_cl, lambda_deriv_cl, alpha_deriv_cl)
        = expressions(calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                      calc_if<!is_constant<T_loc_cl>::value>(mu_deriv),
                      calc_if<!is_constant<T_scale_cl>::value>(lambda_deriv),
                      calc_if<!is_constant<T_shape_cl>::value>(alpha_deriv));

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
  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
#endif
