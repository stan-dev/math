#ifndef STAN_MATH_OPENCL_PRIM_PARETO_CDF_HPP
#define STAN_MATH_OPENCL_PRIM_PARETO_CDF_HPP
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
return_type_t<T_y_cl, T_scale_cl, T_shape_cl> pareto_cdf(
    const T_y_cl& y, const T_scale_cl& y_min, const T_shape_cl& alpha) {
  static const char* function = "pareto_cdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_scale_cl, T_shape_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         y_min, "Scale parameter", alpha);
  const size_t N = max_size(y, y_min, alpha);
  if (N == 0) {
    return 1.0;
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
  auto cond = y_val == INFTY;
  auto log_dbl = log(elt_divide(y_min_val, y_val));
  auto y_min_inv = elt_divide(1.0, y_min_val);
  auto cdf_n = 1.0 - exp(elt_multiply(alpha_val, log_dbl));
  auto cdf_expr = colwise_prod(select(cond, 1.0, cdf_n));
  auto y_y_min_deriv = elt_multiply(alpha_val, y_min_inv);
  auto y_min_alpha_deriv
      = elt_divide(exp(elt_multiply(alpha_val, log_dbl)), cdf_n);

  auto y_deriv = elt_divide(
      elt_multiply(y_y_min_deriv, exp(elt_multiply(alpha_val + 1.0, log_dbl))),
      cdf_n);
  auto y_min_deriv = elt_multiply(-y_y_min_deriv, y_min_alpha_deriv);
  auto alpha_deriv = elt_multiply(y_min_alpha_deriv, -log_dbl);

  matrix_cl<double> cdf_cl;
  matrix_cl<char> any_y_lower_than_y_min_cl;
  matrix_cl<double> y_min_deriv_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> alpha_deriv_cl;

  results(check_y_nonnegative, check_y_min_positive_finite,
          check_alpha_positive_finite, any_y_lower_than_y_min_cl, cdf_cl,
          y_deriv_cl, y_min_deriv_cl, alpha_deriv_cl)
      = expressions(y_not_nonnegative_expr, y_min_positive_finite_expr,
                    alpha_positive_finite_expr, any_y_lower_than_y_min,
                    cdf_expr, calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_scale_cl>::value>(y_min_deriv),
                    calc_if<!is_constant<T_shape_cl>::value>(alpha_deriv));

  if (from_matrix_cl(any_y_lower_than_y_min_cl).maxCoeff()) {
    return 0;
  }

  T_partials_return cdf = from_matrix_cl(cdf_cl).prod();

  operands_and_partials<decltype(y_col), decltype(y_min_col),
                        decltype(alpha_col)>
      ops_partials(y_col, y_min_col, alpha_col);
  if (!is_constant_all<T_y_cl, T_scale_cl, T_shape_cl>::value) {
    results(y_min_deriv_cl, y_deriv_cl, alpha_deriv_cl) = expressions(
        calc_if<!is_constant<T_scale_cl>::value>(y_min_deriv_cl * cdf),
        calc_if<!is_constant<T_y_cl>::value>(y_deriv_cl * cdf),
        calc_if<!is_constant<T_shape_cl>::value>(alpha_deriv_cl * cdf));

    if (!is_constant<T_y_cl>::value) {
      ops_partials.edge1_.partials_ = std::move(y_deriv_cl);
    }
    if (!is_constant<T_scale_cl>::value) {
      ops_partials.edge2_.partials_ = std::move(y_min_deriv_cl);
    }
    if (!is_constant<T_shape_cl>::value) {
      ops_partials.edge3_.partials_ = std::move(alpha_deriv_cl);
    }
  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
#endif
