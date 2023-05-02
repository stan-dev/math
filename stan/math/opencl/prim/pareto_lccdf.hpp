#ifndef STAN_MATH_OPENCL_PRIM_PARETO_LCCDF_HPP
#define STAN_MATH_OPENCL_PRIM_PARETO_LCCDF_HPP
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
 * Returns the Pareto cumulative density function. Given
 * containers of matching sizes, returns the product of probabilities.
 *
 * @tparam T_y_cl type of scalar outcome
 * @tparam T_scale_cl type of location
 * @tparam T_shape_cl type of scale
 * @param y (Sequence of) scalar(s).
 * @param y_min (Sequence of) location(s).
 * @param alpha (Sequence of) scale(s).
 * @return The log of the product of densities.
 */
template <
    typename T_y_cl, typename T_scale_cl, typename T_shape_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_scale_cl,
                                                T_shape_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_scale_cl, T_shape_cl>* = nullptr>
return_type_t<T_y_cl, T_scale_cl, T_shape_cl> pareto_lccdf(
    const T_y_cl& y, const T_scale_cl& y_min, const T_shape_cl& alpha) {
  static const char* function = "pareto_lccdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_scale_cl, T_shape_cl>;
  using std::isfinite;
  using std::isinf;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         y_min, "Scale parameter", alpha);
  const size_t N = max_size(y, y_min, alpha);
  if (N == 0) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& y_min_col = as_column_vector_or_scalar(y_min);
  const auto& alpha_col = as_column_vector_or_scalar(alpha);

  const auto& y_val = value_of(y_col);
  const auto& y_min_val = value_of(y_min_col);
  const auto& alpha_val = value_of(alpha_col);

  auto check_y_nonnegative
      = check_cl(function, "Random variable", y_val, "nonnegative");
  auto y_not_nonnegative_expr = 0 <= y_val;
  auto check_y_min_positive_finite
      = check_cl(function, "Scale parameter", y_min_val, "positive finite");
  auto y_min_positive_finite_expr = 0 < y_min_val && isfinite(y_min_val);
  auto check_alpha_positive_finite
      = check_cl(function, "Shape parameter", alpha_val, "positive finite");
  auto alpha_positive_finite_expr = 0 < alpha_val && isfinite(alpha_val);

  auto any_y_lower_than_y_min = colwise_max(cast<char>(y_val < y_min_val));
  auto any_y_inf = colwise_max(cast<char>(isinf(y_val)));

  auto log_quot = log(elt_divide(y_min_val, y_val));
  auto lccdf_expr = colwise_sum(elt_multiply(alpha_val, log_quot));

  auto alpha_deriv = log_quot;
  auto y_min_deriv = elt_divide(alpha_val, y_min_val);
  auto y_deriv = elt_multiply(-y_min_deriv, exp(log_quot));

  matrix_cl<double> lccdf_cl;
  matrix_cl<char> any_y_lower_than_y_min_cl;
  matrix_cl<char> any_y_inf_cl;
  matrix_cl<double> y_min_deriv_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> alpha_deriv_cl;

  results(check_y_nonnegative, check_y_min_positive_finite,
          check_alpha_positive_finite, any_y_lower_than_y_min_cl, any_y_inf_cl,
          lccdf_cl, y_deriv_cl, y_min_deriv_cl, alpha_deriv_cl)
      = expressions(y_not_nonnegative_expr, y_min_positive_finite_expr,
                    alpha_positive_finite_expr, any_y_lower_than_y_min,
                    any_y_inf, lccdf_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_scale_cl>::value>(y_min_deriv),
                    calc_if<!is_constant<T_shape_cl>::value>(alpha_deriv));

  if (from_matrix_cl(any_y_lower_than_y_min_cl).maxCoeff()) {
    return 0;
  }

  if (from_matrix_cl(any_y_inf_cl).maxCoeff()) {
    return NEGATIVE_INFTY;
  }

  T_partials_return lccdf = from_matrix_cl(lccdf_cl).sum();

  auto ops_partials = make_partials_propagator(y_col, y_min_col, alpha_col);

  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  if (!is_constant<T_scale_cl>::value) {
    partials<1>(ops_partials) = std::move(y_min_deriv_cl);
  }
  if (!is_constant<T_shape_cl>::value) {
    partials<2>(ops_partials) = std::move(alpha_deriv_cl);
  }
  return ops_partials.build(lccdf);
}

}  // namespace math
}  // namespace stan
#endif
#endif
