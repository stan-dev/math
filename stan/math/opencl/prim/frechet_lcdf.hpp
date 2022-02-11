#ifndef STAN_MATH_OPENCL_PRIM_FRECHET_LCDF_HPP
#define STAN_MATH_OPENCL_PRIM_FRECHET_LCDF_HPP
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
 * Returns the frechet log cumulative distribution function for
 * the given location, and scale. If given containers of matching sizes returns
 * the product of probabilities.
 *
 * @tparam T_y_cl type of scalar outcome
 * @tparam T_shape_cl type of location
 * @tparam T_scale_cl type of scale
 * @param y (Sequence of) scalar(s).
 * @param alpha (Sequence of) location(s).
 * @param sigma (Sequence of) scale(s).
 * @return The log of the product of densities.
 */
template <
    typename T_y_cl, typename T_shape_cl, typename T_scale_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_shape_cl,
                                                T_scale_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_shape_cl, T_scale_cl>* = nullptr>
return_type_t<T_y_cl, T_shape_cl, T_scale_cl> frechet_lcdf(
    const T_y_cl& y, const T_shape_cl& alpha, const T_scale_cl& sigma) {
  static const char* function = "frechet_lcdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_shape_cl, T_scale_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Shape parameter",
                         alpha, "Scale parameter", sigma);
  const size_t N = max_size(y, alpha, sigma);
  if (N == 0) {
    return 1.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& alpha_col = as_column_vector_or_scalar(alpha);
  const auto& sigma_col = as_column_vector_or_scalar(sigma);

  const auto& y_val = value_of(y_col);
  const auto& alpha_val = value_of(alpha_col);
  const auto& sigma_val = value_of(sigma_col);

  auto check_y_positive
      = check_cl(function, "Random variable", y_val, "positive");
  auto y_positive = y_val > 0;
  auto check_alpha_positive_finite
      = check_cl(function, "Shape parameter", alpha_val, "positive finite");
  auto alpha_positive_finite_expr = alpha_val > 0 && isfinite(alpha_val);
  auto check_sigma_positive_finite
      = check_cl(function, "Scale parameter", sigma_val, "positive finite");
  auto sigma_positive_finite_expr = 0 < sigma_val && isfinite(sigma_val);

  auto pow_n = pow(elt_divide(sigma_val, y_val), alpha_val);
  auto lcdf_expr = colwise_sum(pow_n);

  auto y_deriv = elt_multiply(elt_divide(alpha_val, y_val), pow_n);
  auto alpha_deriv = elt_multiply(log(elt_divide(y_val, sigma_val)), pow_n);
  auto sigma_deriv = elt_multiply(elt_divide(alpha_val, -sigma_val), pow_n);

  matrix_cl<double> lcdf_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> alpha_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;

  results(check_y_positive, check_alpha_positive_finite,
          check_sigma_positive_finite, lcdf_cl, y_deriv_cl, alpha_deriv_cl,
          sigma_deriv_cl)
      = expressions(y_positive, alpha_positive_finite_expr,
                    sigma_positive_finite_expr, lcdf_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_shape_cl>::value>(alpha_deriv),
                    calc_if<!is_constant<T_scale_cl>::value>(sigma_deriv));

  T_partials_return lcdf = -(from_matrix_cl(lcdf_cl)).sum();

  operands_and_partials<decltype(y_col), decltype(alpha_col),
                        decltype(sigma_col)>
      ops_partials(y_col, alpha_col, sigma_col);

  if (!is_constant<T_y_cl>::value) {
    ops_partials.edge1_.partials_ = std::move(y_deriv_cl);
  }
  if (!is_constant<T_shape_cl>::value) {
    ops_partials.edge2_.partials_ = std::move(alpha_deriv_cl);
  }
  if (!is_constant<T_scale_cl>::value) {
    ops_partials.edge3_.partials_ = std::move(sigma_deriv_cl);
  }
  return ops_partials.build(lcdf);
}

}  // namespace math
}  // namespace stan
#endif
#endif
