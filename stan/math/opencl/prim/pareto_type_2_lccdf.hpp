#ifndef STAN_MATH_OPENCL_PRIM_DOUBLE_PARETO_TYPE_2_LCCDF_HPP
#define STAN_MATH_OPENCL_PRIM_DOUBLE_PARETO_TYPE_2_LCCDF_HPP
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
 * Returns the pareto type 2 log complementaty cumulative density function.
 * Given containers of matching sizes, returns the sum of log probabilities.
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
return_type_t<T_y_cl, T_loc_cl, T_scale_cl, T_shape_cl> pareto_type_2_lccdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_scale_cl& lambda,
    const T_shape_cl& alpha) {
  static const char* function = "pareto_type_2_lccdf(OpenCL)";
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

  auto log_temp = log1p(elt_divide(diff, lambda_val));
  auto lccdf_expr = colwise_sum(elt_multiply(alpha_val, log_temp));

  auto mu_deriv = elt_divide(alpha_val, diff + lambda_val);
  auto y_deriv = -mu_deriv;
  auto lambda_deriv = elt_divide(elt_multiply(mu_deriv, diff), lambda_val);
  auto alpha_deriv = -log_temp;

  matrix_cl<double> lccdf_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> lambda_deriv_cl;
  matrix_cl<double> alpha_deriv_cl;

  results(check_y_nonnegative, check_lambda_positive_finite,
          check_alpha_positive_finite, check_diff_nonnegative, lccdf_cl,
          y_deriv_cl, mu_deriv_cl, lambda_deriv_cl, alpha_deriv_cl)
      = expressions(y_nonnegative_expr, lambda_positive_finite_expr,
                    alpha_positive_finite_expr, diff_nonnegative_expr,
                    lccdf_expr, calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_loc_cl>::value>(mu_deriv),
                    calc_if<!is_constant<T_scale_cl>::value>(lambda_deriv),
                    calc_if<!is_constant<T_shape_cl>::value>(alpha_deriv));

  T_partials_return lccdf = -from_matrix_cl(lccdf_cl).sum();

  operands_and_partials<decltype(y_col), decltype(mu_col), decltype(lambda_col),
                        decltype(alpha_col)>
      ops_partials(y_col, mu_col, lambda_col, alpha_col);
  if (!is_constant_all<T_y_cl, T_loc_cl, T_scale_cl, T_shape_cl>::value) {
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
  return ops_partials.build(lccdf);
}

}  // namespace math
}  // namespace stan
#endif
#endif
