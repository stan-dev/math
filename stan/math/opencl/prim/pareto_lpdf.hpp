#ifndef STAN_MATH_OPENCL_PRIM_PARETO_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_PARETO_LPDF_HPP
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
 * The log of the Cauchy density for the specified scalar(s) given
 * the specified location parameter(s) and scale parameter(s). y,
 * y_min, or alpha can each either be scalar a vector.  Any vector
 * inputs must be the same length.
 *
 * <p> The result log probability is defined to be the sum of
 * the log probabilities for each observation/y_min/alpha triple.
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
    bool propto, typename T_y_cl, typename T_scale_cl, typename T_shape_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_scale_cl,
                                                T_shape_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_scale_cl, T_shape_cl>* = nullptr>
return_type_t<T_y_cl, T_scale_cl, T_shape_cl> pareto_lpdf(
    const T_y_cl& y, const T_scale_cl& y_min, const T_shape_cl& alpha) {
  static constexpr const char* function = "pareto_lpdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_scale_cl, T_shape_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Scale parameter",
                         y_min, "Shape parameter", alpha);
  const size_t N = max_size(y, y_min, alpha);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_scale_cl, T_shape_cl>::value) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& y_min_col = as_column_vector_or_scalar(y_min);
  const auto& alpha_col = as_column_vector_or_scalar(alpha);

  const auto& y_val = value_of(y_col);
  const auto& y_min_val = value_of(y_min_col);
  const auto& alpha_val = value_of(alpha_col);

  auto check_y_not_nan
      = check_cl(function, "Random variable", y_val, "not NaN");
  auto y_not_nan = !isnan(y_val);
  auto check_y_min_positive_finite
      = check_cl(function, "Scale parameter", y_min_val, "positive finite");
  auto y_min_positive_finite = 0 < y_min_val && isfinite(y_min_val);
  auto check_alpha_positive_finite
      = check_cl(function, "Shape parameter", alpha_val, "positive finite");
  auto alpha_positive_finite = 0 < alpha_val && isfinite(alpha_val);

  auto y_less_than_y_min = colwise_max(cast<char>(y_val < y_min_val));
  auto log_y = log(y_val);
  auto inv_y = elt_divide(1.0, y_val);
  auto log_y_min = log(y_min_val);
  auto logp1 = static_select<include_summand<propto, T_shape_cl>::value>(
      log(alpha_val), constant(0.0, N, 1));
  auto logp2
      = static_select<include_summand<propto, T_y_cl, T_shape_cl>::value>(
          logp1 - elt_multiply(alpha_val, log_y) - log_y, logp1);
  auto logp_expr = colwise_sum(
      static_select<include_summand<propto, T_scale_cl, T_shape_cl>::value>(
          logp2 + elt_multiply(alpha_val, log_y_min), logp2));

  auto y_deriv = -(elt_multiply(alpha_val, inv_y) + inv_y);
  auto y_min_deriv = elt_divide(alpha_val, y_min_val);
  auto alpha_deriv = elt_divide(1.0, alpha_val) + log_y_min - log_y;

  matrix_cl<char> y_less_than_y_min_cl;
  matrix_cl<double> logp_cl;
  matrix_cl<double> y_min_deriv_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> alpha_deriv_cl;

  results(check_y_not_nan, check_y_min_positive_finite,
          check_alpha_positive_finite, y_less_than_y_min_cl, logp_cl,
          y_deriv_cl, y_min_deriv_cl, alpha_deriv_cl)
      = expressions(y_not_nan, y_min_positive_finite, alpha_positive_finite,
                    y_less_than_y_min, logp_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_scale_cl>::value>(y_min_deriv),
                    calc_if<!is_constant<T_shape_cl>::value>(alpha_deriv));

  if (from_matrix_cl(y_less_than_y_min_cl).any()) {
    return LOG_ZERO;
  }

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

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
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
