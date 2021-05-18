#ifndef STAN_MATH_OPENCL_PRIM_INV_GAMMA_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_INV_GAMMA_LPDF_HPP
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
 * The log of an inverse gamma density for y with the specified
 * shape and scale parameters.
 * Shape and scale parameters must be greater than 0.
 * y must be greater than 0.
 *
 * @param y A scalar variable.
 * @param alpha Shape parameter.
 * @param beta Scale parameter.
 * @throw std::domain_error if alpha is not greater than 0.
 * @throw std::domain_error if beta is not greater than 0.
 * @throw std::domain_error if y is not greater than 0.
 * @tparam T_y Type of scalar.
 * @tparam T_shape Type of shape.
 * @tparam T_scale Type of scale.
 */
template <
    bool propto, typename T_y_cl, typename T_shape_cl, typename T_scale_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_shape_cl,
                                                T_scale_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_shape_cl, T_scale_cl>* = nullptr>
return_type_t<T_y_cl, T_shape_cl, T_scale_cl> inv_gamma_lpdf(
    const T_y_cl& y, const T_shape_cl& alpha, const T_scale_cl& beta) {
  using std::isfinite;
  using std::isnan;
  static const char* function = "inv_gamma_lpdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_shape_cl, T_scale_cl>;

  check_consistent_sizes(function, "Random variable", y,
                         "First shape parameter", alpha,
                         "Second shape parameter", beta);
  const size_t N = max_size(y, alpha, beta);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_shape_cl, T_scale_cl>::value) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& alpha_col = as_column_vector_or_scalar(alpha);
  const auto& beta_col = as_column_vector_or_scalar(beta);

  const auto& y_val = value_of(y_col);
  const auto& alpha_val = value_of(alpha_col);
  const auto& beta_val = value_of(beta_col);

  operands_and_partials<decltype(y_col), decltype(alpha_col),
                        decltype(beta_col)>
      ops_partials(y_col, alpha_col, beta_col);

  auto check_y_not_nan
      = check_cl(function, "Random variable", y_val, "not NaN");
  auto y_not_nan = !isnan(y_val);
  auto check_alpha_pos_finite
      = check_cl(function, "Shape parameter", alpha_val, "positive finite");
  auto alpha_pos_finite = alpha_val > 0 && isfinite(alpha_val);
  auto check_beta_pos_finite
      = check_cl(function, "Scale parameter", beta_val, "positive finite");
  auto beta_pos_finite = beta_val > 0 && isfinite(beta_val);

  auto any_y_nonpositive = colwise_max(cast<char>(y_val <= 0));
  auto log_y = log(y_val);
  auto log_beta = log(beta_val);
  auto inv_y = elt_divide(1.0, y_val);

  auto logp1 = static_select<include_summand<propto, T_shape_cl>::value>(
      -lgamma(alpha_val), constant(0.0, N, 1));
  auto logp2
      = static_select<include_summand<propto, T_shape_cl, T_scale_cl>::value>(
          logp1 + elt_multiply(alpha_val, log_beta), logp1);
  auto logp3
      = static_select<include_summand<propto, T_y_cl, T_shape_cl>::value>(
          logp2 - elt_multiply(alpha_val + 1.0, log_y), logp2);
  auto logp_expr = colwise_sum(
      static_select<include_summand<propto, T_y_cl, T_scale_cl>::value>(
          logp3 - elt_multiply(beta_val, inv_y), logp3));

  auto y_deriv
      = elt_multiply(elt_multiply(beta_val, inv_y) - alpha_val - 1, inv_y);
  auto alpha_deriv = log_beta - digamma(alpha_val) - log_y;
  auto beta_deriv = elt_divide(alpha_val, beta_val) - inv_y;

  matrix_cl<char> any_y_nonpositive_cl;
  matrix_cl<double> logp_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> alpha_deriv_cl;
  matrix_cl<double> beta_deriv_cl;

  results(check_alpha_pos_finite, check_beta_pos_finite, check_y_not_nan,
          any_y_nonpositive_cl, logp_cl, y_deriv_cl, alpha_deriv_cl,
          beta_deriv_cl)
      = expressions(alpha_pos_finite, beta_pos_finite, y_not_nan,
                    any_y_nonpositive, logp_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_shape_cl>::value>(alpha_deriv),
                    calc_if<!is_constant<T_scale_cl>::value>(beta_deriv));

  if (from_matrix_cl(any_y_nonpositive_cl).any()) {
    return LOG_ZERO;
  }

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  if (!is_constant<T_y_cl>::value) {
    ops_partials.edge1_.partials_ = std::move(y_deriv_cl);
  }
  if (!is_constant<T_shape_cl>::value) {
    ops_partials.edge2_.partials_ = std::move(alpha_deriv_cl);
  }
  if (!is_constant<T_scale_cl>::value) {
    ops_partials.edge3_.partials_ = std::move(beta_deriv_cl);
  }

  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
