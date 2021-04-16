#ifndef STAN_MATH_OPENCL_PRIM_UNIFORM_LCCDF_HPP
#define STAN_MATH_OPENCL_PRIM_UNIFORM_LCCDF_HPP
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
 * Returns the log uniform complementary cumulative distribution function for
 * the given location, and scale. If given containers of matching sizes returns
 * the product of probabilities.
 *
 * @tparam T_y_cl type of scalar outcome
 * @tparam T_low_cl type of location
 * @tparam T_high_cl type of scale
 * @param y (Sequence of) scalar(s).
 * @param alpha (Sequence of) location(s).
 * @param beta (Sequence of) scale(s).
 * @return The log of the product of densities.
 */
template <typename T_y_cl, typename T_low_cl, typename T_high_cl,
          require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_low_cl,
                                                      T_high_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl, T_low_cl, T_high_cl>* = nullptr>
return_type_t<T_y_cl, T_low_cl, T_high_cl> uniform_lccdf(
    const T_y_cl& y, const T_low_cl& alpha, const T_high_cl& beta) {
  static const char* function = "uniform_lccdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_low_cl, T_high_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         alpha, "Scale parameter", beta);
  const size_t N = max_size(y, alpha, beta);
  if (N == 0) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& alpha_col = as_column_vector_or_scalar(alpha);
  const auto& beta_col = as_column_vector_or_scalar(beta);

  const auto& y_val = value_of(y_col);
  const auto& alpha_val = value_of(alpha_col);
  const auto& beta_val = value_of(beta_col);

  auto check_y_not_nan
      = check_cl(function, "Random variable", y_val, "not NaN");
  auto y_not_nan_expr = !isnan(y_val);
  auto check_alpha_finite
      = check_cl(function, "Lower bound parameter", alpha_val, "finite");
  auto alpha_finite_expr = isfinite(alpha_val);
  auto check_beta_finite
      = check_cl(function, "Upper bound parameter", beta_val, "finite");
  auto beta_finite_expr = isfinite(beta_val);
  auto b_minus_a = beta_val - alpha_val;
  auto check_diff_positive = check_cl(
      function, "Difference between upper and lower bound parameters", beta_val,
      "positive");
  auto diff_positive_expr = b_minus_a > 0.0;

  auto any_y_out_of_bounds
      = colwise_max(cast<char>(y_val < alpha_val || y_val > beta_val));
  auto y_minus_alpha = y_val - alpha_val;
  auto ccdf_n = 1.0 - elt_divide(y_minus_alpha, b_minus_a);
  auto lccdf_expr = colwise_sum(log(ccdf_n));

  auto y_deriv = elt_divide(1.0, -elt_multiply(b_minus_a, ccdf_n));
  auto rep_deriv = elt_divide(1.0, elt_multiply(square(b_minus_a), ccdf_n));
  auto low_deriv = elt_multiply(beta_val - y_val, rep_deriv);
  auto high_deriv = elt_multiply(y_minus_alpha, rep_deriv);

  matrix_cl<char> any_y_out_of_bounds_cl;
  matrix_cl<double> lccdf_cl;
  matrix_cl<double> alpha_deriv_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> beta_deriv_cl;

  results(check_y_not_nan, check_alpha_finite, check_beta_finite,
          check_diff_positive, any_y_out_of_bounds_cl, lccdf_cl, y_deriv_cl,
          alpha_deriv_cl, beta_deriv_cl)
      = expressions(y_not_nan_expr, alpha_finite_expr, beta_finite_expr,
                    diff_positive_expr, any_y_out_of_bounds, lccdf_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_low_cl>::value>(low_deriv),
                    calc_if<!is_constant<T_high_cl>::value>(high_deriv));

  if (from_matrix_cl(any_y_out_of_bounds_cl).maxCoeff()) {
    return 0.0;
  }

  T_partials_return lccdf = from_matrix_cl(lccdf_cl).sum();

  operands_and_partials<decltype(y_col), decltype(alpha_col),
                        decltype(beta_col)>
      ops_partials(y_col, alpha_col, beta_col);

  if (!is_constant<T_y_cl>::value) {
    ops_partials.edge1_.partials_ = std::move(y_deriv_cl);
  }
  if (!is_constant<T_low_cl>::value) {
    ops_partials.edge2_.partials_ = std::move(alpha_deriv_cl);
  }
  if (!is_constant<T_high_cl>::value) {
    ops_partials.edge3_.partials_ = std::move(beta_deriv_cl);
  }
  return ops_partials.build(lccdf);
}

}  // namespace math
}  // namespace stan
#endif
#endif
