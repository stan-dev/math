#ifndef STAN_MATH_OPENCL_PRIM_GAMMA_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_GAMMA_LPDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * The log of a gamma density for y with the specified
 * shape and inverse scale parameters.
 * Shape and inverse scale parameters must be greater than 0.
 * y must be greater than or equal to 0.
 *
 \f{eqnarray*}{
 y &\sim& \mbox{\sf{Gamma}}(\alpha, \beta) \\
 \log (p (y \, |\, \alpha, \beta) ) &=& \log \left(
 \frac{\beta^\alpha}{\Gamma(\alpha)} y^{\alpha - 1} \exp^{- \beta y} \right) \\
 &=& \alpha \log(\beta) - \log(\Gamma(\alpha)) + (\alpha - 1) \log(y) - \beta
 y\\ & & \mathrm{where} \; y > 0 \f}
 *
 * @tparam T_y_cl type of scalar
 * @tparam T_shape_cl type of shape
 * @tparam T_inv_scale_cl type of inverse scale
 * @param y A scalar variable.
 * @param alpha Shape parameter.
 * @param beta Inverse scale parameter.
 * @throw std::domain_error if alpha is not greater than 0.
 * @throw std::domain_error if beta is not greater than 0.
 * @throw std::domain_error if y is not greater than or equal to 0.
 */
template <bool propto, typename T_y_cl, typename T_shape_cl,
          typename T_inv_scale_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_y_cl, T_shape_cl, T_inv_scale_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl, T_shape_cl,
                                        T_inv_scale_cl>* = nullptr>
return_type_t<T_y_cl, T_shape_cl, T_inv_scale_cl> gamma_lpdf(
    const T_y_cl& y, const T_shape_cl& alpha, const T_inv_scale_cl& beta) {
  using std::isfinite;
  using std::isnan;
  static const char* function = "gamma_lpdf(OpenCL)";
  using T_partials_return
      = partials_return_t<T_y_cl, T_shape_cl, T_inv_scale_cl>;

  check_consistent_sizes(function, "Random variable", y, "Shape parameter",
                         alpha, "Inverse scale parameter", beta);
  const size_t N = max_size(y, alpha, beta);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_shape_cl, T_inv_scale_cl>::value) {
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
  auto y_not_nan_expr = y_val > 0 && isfinite(y_val);
  auto check_alpha_pos_finite
      = check_cl(function, "Shape parameter", alpha_val, "positive finite");
  auto alpha_pos_finite_expr = alpha_val > 0 && isfinite(alpha_val);
  auto check_beta_pos_finite = check_cl(function, "Inverse scale parameter",
                                        beta_val, "positive finite");
  auto beta_pos_finite_expr = beta_val > 0 && isfinite(beta_val);

  auto any_y_negative_expr = colwise_max(cast<char>(y_val < 0));
  auto log_y_expr = log(y_val);
  auto log_beta_expr = log(beta_val);
  auto logp1_expr = static_select<include_summand<propto, T_shape_cl>::value>(
      -lgamma(alpha_val), constant(0.0, N, 1));
  auto logp2_expr = static_select<
      include_summand<propto, T_shape_cl, T_inv_scale_cl>::value>(
      logp1_expr + elt_multiply(alpha_val, log_beta_expr), logp1_expr);
  auto logp3_expr
      = static_select<include_summand<propto, T_y_cl, T_shape_cl>::value>(
          logp2_expr + elt_multiply(alpha_val - 1.0, log_y_expr), logp2_expr);
  auto logp_expr = colwise_sum(
      static_select<include_summand<propto, T_y_cl, T_inv_scale_cl>::value>(
          logp3_expr - elt_multiply(beta_val, y_val), logp3_expr));

  auto y_deriv_expr = elt_divide(alpha_val - 1, y_val) - beta_val;
  auto alpha_deriv_expr = log_beta_expr + log_y_expr - digamma(alpha_val);
  auto beta_deriv_expr = elt_divide(alpha_val, beta_val) - y_val;

  matrix_cl<char> any_y_negative_cl;
  matrix_cl<double> logp_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> alpha_deriv_cl;
  matrix_cl<double> beta_deriv_cl;

  results(check_y_not_nan, check_alpha_pos_finite, check_beta_pos_finite,
          any_y_negative_cl, logp_cl, y_deriv_cl, alpha_deriv_cl, beta_deriv_cl)
      = expressions(
          y_not_nan_expr, alpha_pos_finite_expr, beta_pos_finite_expr,
          any_y_negative_expr, logp_expr,
          calc_if<!is_constant<T_y_cl>::value>(y_deriv_expr),
          calc_if<!is_constant<T_shape_cl>::value>(alpha_deriv_expr),
          calc_if<!is_constant<T_inv_scale_cl>::value>(beta_deriv_expr));

  if (from_matrix_cl(any_y_negative_cl).any()) {
    return LOG_ZERO;
  }

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  operands_and_partials<decltype(y_col), decltype(alpha_col),
                        decltype(beta_col)>
      ops_partials(y_col, alpha_col, beta_col);
  if (!is_constant<T_y_cl>::value) {
    ops_partials.edge1_.partials_ = std::move(y_deriv_cl);
  }
  if (!is_constant<T_shape_cl>::value) {
    ops_partials.edge2_.partials_ = std::move(alpha_deriv_cl);
  }
  if (!is_constant<T_inv_scale_cl>::value) {
    ops_partials.edge3_.partials_ = std::move(beta_deriv_cl);
  }

  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
