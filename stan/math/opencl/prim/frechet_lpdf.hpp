#ifndef STAN_MATH_OPENCL_PRIM_FRECHET_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_FRECHET_LPDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * The log of the frechet density for the specified scalar(s) given the
 * specified sample stan::math::size(s). y, alpha, or sigma can each either be
 * scalar or a vector on OpenCL device. Any vector inputs must be the same
 * length.
 *
 * <p> The result log probability is defined to be the sum of
 * the log probabilities for each observation/alpha/sigma triple.
 *
 * @tparam T_y_cl type of scalar outcome
 * @tparam T_shape_cl type of shape
 * @tparam T_scale_cl type of scale
 * @param y (Sequence of) scalar(s).
 * @param alpha shape
 * @param sigma scale
 * @return The log of the product of densities.
 */
template <
    bool propto, typename T_y_cl, typename T_shape_cl, typename T_scale_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_shape_cl,
                                                T_scale_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_shape_cl, T_scale_cl>* = nullptr>
return_type_t<T_y_cl, T_shape_cl, T_scale_cl> frechet_lpdf(
    const T_y_cl& y, const T_shape_cl& alpha, const T_scale_cl& sigma) {
  using std::isfinite;
  static const char* function = "frechet_lpdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_shape_cl, T_scale_cl>;

  check_consistent_sizes(function, "Random variable", y, "Shape parameter",
                         alpha, "Scale parameter", sigma);
  const size_t N = max_size(y, alpha, sigma);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_shape_cl, T_scale_cl>::value) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& alpha_col = as_column_vector_or_scalar(alpha);
  const auto& sigma_col = as_column_vector_or_scalar(sigma);

  const auto& y_val = value_of(y_col);
  const auto& alpha_val = value_of(alpha_col);
  const auto& sigma_val = value_of(sigma_col);

  auto ops_partials = make_partials_propagator(y_col, alpha_col, sigma_col);

  auto check_y_positive
      = check_cl(function, "Random variable", y_val, "positive");
  auto y_positive_expr = 0 < y_val;
  auto check_alpha_pos_finite
      = check_cl(function, "Shape parameter", alpha_val, "positive finite");
  auto alpha_pos_finite_expr = alpha_val > 0 && isfinite(alpha_val);
  auto check_sigma_pos_finite
      = check_cl(function, "Scale parameter", sigma_val, "positive finite");
  auto sigma_pos_finite_expr = sigma_val > 0 && isfinite(sigma_val);

  auto log_y_expr = log(y_val);
  auto sigma_div_y_pow_alpha_expr
      = pow(elt_divide(sigma_val, y_val), alpha_val);
  auto log_sigma_expr = log(sigma_val);
  auto logp1_expr = -sigma_div_y_pow_alpha_expr;
  auto logp2_expr = static_select<include_summand<propto, T_shape_cl>::value>(
      logp1_expr + log(alpha_val), logp1_expr);
  auto logp3_expr
      = static_select<include_summand<propto, T_y_cl, T_shape_cl>::value>(
          logp2_expr - elt_multiply(alpha_val + 1.0, log_y_expr), logp2_expr);
  auto logp_expr = colwise_sum(
      static_select<include_summand<propto, T_shape_cl, T_scale_cl>::value>(
          logp3_expr + elt_multiply(alpha_val, log_sigma_expr), logp3_expr));

  auto alpha_deriv_expr = elt_divide(1.0, alpha_val)
                          + elt_multiply(1 - sigma_div_y_pow_alpha_expr,
                                         log_sigma_expr - log_y_expr);
  auto y_deriv_expr = elt_divide(
      elt_multiply(alpha_val, sigma_div_y_pow_alpha_expr) - alpha_val - 1,
      y_val);
  auto sigma_deriv_expr = elt_multiply(elt_divide(alpha_val, sigma_val),
                                       1 - sigma_div_y_pow_alpha_expr);

  matrix_cl<double> logp_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> alpha_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;

  results(check_y_positive, check_alpha_pos_finite, check_sigma_pos_finite,
          logp_cl, y_deriv_cl, alpha_deriv_cl, sigma_deriv_cl)
      = expressions(y_positive_expr, alpha_pos_finite_expr,
                    sigma_pos_finite_expr, logp_expr, y_deriv_expr,
                    alpha_deriv_expr, sigma_deriv_expr);

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  if (!is_constant<T_shape_cl>::value) {
    partials<1>(ops_partials) = std::move(alpha_deriv_cl);
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
