#ifndef STAN_MATH_OPENCL_PRIM_LOGISTIC_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_LOGISTIC_LPDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <stan/math/prim/fun/exp.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * The log of a logistic density for y with the specified
 * location and scale parameters.
 * Shape and scale parameters must be greater than 0.
 * y must be greater than 0.
 *
 * @param y A scalar variable.
 * @param mu Shape parameter.
 * @param sigma Scale parameter.
 * @throw std::domain_error if mu is not greater than 0.
 * @throw std::domain_error if sigma is not greater than 0.
 * @throw std::domain_error if y is not greater than 0.
 * @tparam T_y_cl Type of scalar.
 * @tparam T_loc_cl Type of shape.
 * @tparam T_scale_cl Type of scale.
 */
template <
    bool propto, typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_loc_cl,
                                                T_scale_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_scale_cl>* = nullptr>
return_type_t<T_y_cl, T_loc_cl, T_scale_cl> logistic_lpdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_scale_cl& sigma) {
  using std::isfinite;
  static const char* function = "logistic_lpdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_loc_cl, T_scale_cl>;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
  const size_t N = max_size(y, mu, sigma);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_loc_cl, T_scale_cl>::value) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& mu_col = as_column_vector_or_scalar(mu);
  const auto& sigma_col = as_column_vector_or_scalar(sigma);

  const auto& y_val = value_of(y_col);
  const auto& mu_val = value_of(mu_col);
  const auto& sigma_val = value_of(sigma_col);

  operands_and_partials<decltype(y_col), decltype(mu_col), decltype(sigma_col)>
      ops_partials(y_col, mu_col, sigma_col);

  auto check_y_finite = check_cl(function, "Random variable", y_val, "finite");
  auto y_finite = isfinite(y_val);
  auto check_mu_finite
      = check_cl(function, "Location parameter", mu_val, "finite");
  auto mu_finite = isfinite(mu_val);
  auto check_sigma_pos_finite
      = check_cl(function, "Scale parameter", sigma_val, "positive finite");
  auto sigma_pos_finite = sigma_val > 0.0 && isfinite(sigma_val);

  auto inv_sigma = elt_divide(1.0, sigma_val);
  auto y_minus_mu = y_val - mu_val;
  auto y_minus_mu_div_sigma = elt_multiply(y_minus_mu, inv_sigma);

  auto logp1 = -y_minus_mu_div_sigma - 2.0 * log1p(exp(-y_minus_mu_div_sigma));
  auto logp_expr
      = colwise_sum(static_select<include_summand<propto, T_scale_cl>::value>(
          logp1 - log(sigma_val), logp1));

  auto y_deriv = elt_multiply(
      elt_divide(2.0, 1.0 + exp(y_minus_mu_div_sigma)) - 1.0, inv_sigma);
  auto exp_mu_div_sigma = exp(elt_multiply(mu_val, inv_sigma));
  auto mu_deriv = elt_multiply(
      1.0
          - 2.0
                * elt_divide(
                      exp_mu_div_sigma,
                      exp_mu_div_sigma + exp(elt_multiply(y_val, inv_sigma))),
      inv_sigma);
  auto sigma_deriv
      = elt_multiply(-elt_multiply(y_deriv, y_minus_mu) - 1.0, inv_sigma);

  matrix_cl<double> logp_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;

  results(check_mu_finite, check_sigma_pos_finite, check_y_finite, logp_cl,
          y_deriv_cl, mu_deriv_cl, sigma_deriv_cl)
      = expressions(mu_finite, sigma_pos_finite, y_finite, logp_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_loc_cl>::value>(mu_deriv),
                    calc_if<!is_constant<T_scale_cl>::value>(sigma_deriv));

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  if (!is_constant<T_y_cl>::value) {
    ops_partials.edge1_.partials_ = std::move(y_deriv_cl);
  }
  if (!is_constant<T_loc_cl>::value) {
    ops_partials.edge2_.partials_ = std::move(mu_deriv_cl);
  }
  if (!is_constant<T_scale_cl>::value) {
    ops_partials.edge3_.partials_ = std::move(sigma_deriv_cl);
  }

  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
