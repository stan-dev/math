#ifndef STAN_MATH_OPENCL_PRIM_UNIFORM_CDF_HPP
#define STAN_MATH_OPENCL_PRIM_UNIFORM_CDF_HPP
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
 * Returns the uniform cumulative distribution function for the given
 * location, and scale. If given containers of matching sizes
 * returns the product of probabilities.
 *
 * @tparam T_y_cl type of scalar outcome
 * @tparam T_low_cl type of low bounds
 * @tparam T_high_cl type of high bounds
 * @param y (Sequence of) scalar(s).
 * @param alpha Sequence of low bounds.
 * @param beta Sequence of high bounds.
 * @return The product of densities.
 */
template <typename T_y_cl, typename T_low_cl, typename T_high_cl,
          require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_low_cl,
                                                      T_high_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl, T_low_cl, T_high_cl>* = nullptr>
return_type_t<T_y_cl, T_low_cl, T_high_cl> uniform_cdf(const T_y_cl& y,
                                                       const T_low_cl& alpha,
                                                       const T_high_cl& beta) {
  static constexpr const char* function = "uniform_cdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_low_cl, T_high_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         alpha, "Scale parameter", beta);
  const size_t N = max_size(y, alpha, beta);
  if (N == 0) {
    return 1.0;
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
  auto cdf_n = elt_divide(y_val - alpha_val, b_minus_a);
  auto cdf_expr = colwise_prod(cdf_n);

  auto high_deriv1 = elt_divide(1.0, b_minus_a);
  auto y_deriv1 = elt_divide(high_deriv1, cdf_n);
  auto low_deriv1
      = elt_multiply(y_val - beta_val, elt_divide(y_deriv1, b_minus_a));

  matrix_cl<char> any_y_out_of_bounds_cl;
  matrix_cl<double> cdf_cl;
  matrix_cl<double> alpha_deriv_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> beta_deriv_cl;

  results(check_y_not_nan, check_alpha_finite, check_beta_finite,
          check_diff_positive, any_y_out_of_bounds_cl, cdf_cl, y_deriv_cl,
          alpha_deriv_cl, beta_deriv_cl)
      = expressions(y_not_nan_expr, alpha_finite_expr, beta_finite_expr,
                    diff_positive_expr, any_y_out_of_bounds, cdf_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv1),
                    calc_if<!is_constant<T_low_cl>::value>(low_deriv1),
                    calc_if<!is_constant<T_high_cl>::value>(high_deriv1));

  if (from_matrix_cl(any_y_out_of_bounds_cl).maxCoeff()) {
    return 0.0;
  }

  T_partials_return cdf = (from_matrix_cl(cdf_cl)).prod();

  auto alpha_deriv = alpha_deriv_cl * cdf;
  auto y_deriv = y_deriv_cl * cdf;
  auto beta_deriv = beta_deriv_cl * -cdf;

  results(alpha_deriv_cl, y_deriv_cl, beta_deriv_cl)
      = expressions(calc_if<!is_constant<T_low_cl>::value>(alpha_deriv),
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_high_cl>::value>(beta_deriv));

  auto ops_partials = make_partials_propagator(y_col, alpha_col, beta_col);

  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  if (!is_constant<T_low_cl>::value) {
    partials<1>(ops_partials) = std::move(alpha_deriv_cl);
  }
  if (!is_constant<T_high_cl>::value) {
    partials<2>(ops_partials) = std::move(beta_deriv_cl);
  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
#endif
