#ifndef STAN_MATH_OPENCL_PRIM_RAYLEIGH_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_RAYLEIGH_LPDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * The log of an Rayleigh density for y with the specified
 * scale parameter.
 * y and scale parameter must be greater than 0.
 *
 * @tparam T_y_cl type of scalar
 * @tparam T_scale_cl type of inverse scale
 * @param y A scalar variable.
 * @param sigma Inverse scale parameter.
 * @throw std::domain_error if sigma is not greater than 0.
 * @throw std::domain_error if y is not greater than or equal to 0.
 */
template <
    bool propto, typename T_y_cl, typename T_scale_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_scale_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_scale_cl>* = nullptr>
return_type_t<T_y_cl, T_scale_cl> rayleigh_lpdf(const T_y_cl& y,
                                                const T_scale_cl& sigma) {
  static const char* function = "rayleigh_lpdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_scale_cl>;

  check_consistent_sizes(function, "Random variable", y, "Scale parameter",
                         sigma);
  const size_t N = max_size(y, sigma);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_scale_cl>::value) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& sigma_col = as_column_vector_or_scalar(sigma);

  const auto& y_val = value_of(y_col);
  const auto& sigma_val = value_of(sigma_col);

  operands_and_partials<decltype(y_col), decltype(sigma_col)> ops_partials(
      y_col, sigma_col);

  auto check_y_positive
      = check_cl(function, "Random variable", y_val, "positive");
  auto y_positive = y_val > 0;
  auto check_sigma_positive
      = check_cl(function, "Scale parameter", sigma_val, "positive");
  auto sigma_positive = sigma_val > 0;

  auto inv_sigma = elt_divide(1.0, sigma_val);
  auto y_over_sigma = elt_divide(y_val, sigma_val);

  auto logp1 = -0.5 * elt_multiply(y_over_sigma, y_over_sigma);
  auto logp2 = static_select<include_summand<propto, T_scale_cl>::value>(
      logp1 - 2.0 * log(sigma_val), logp1);
  auto logp_expr
      = colwise_sum(static_select<include_summand<propto, T_y_cl>::value>(
          logp2 + log(y_val), logp2));

  auto scaled_diff = elt_multiply(inv_sigma, y_over_sigma);
  auto y_deriv_expr = elt_divide(1.0, y_val) - scaled_diff;
  auto sigma_deriv_expr
      = elt_multiply(y_over_sigma, scaled_diff) - 2.0 * inv_sigma;

  matrix_cl<double> logp_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;

  results(check_y_positive, check_sigma_positive, logp_cl, y_deriv_cl,
          sigma_deriv_cl)
      = expressions(y_positive, sigma_positive, logp_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv_expr),
                    calc_if<!is_constant<T_scale_cl>::value>(sigma_deriv_expr));

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  if (!is_constant<T_y_cl>::value) {
    ops_partials.edge1_.partials_ = std::move(y_deriv_cl);
  }
  if (!is_constant<T_scale_cl>::value) {
    ops_partials.edge2_.partials_ = std::move(sigma_deriv_cl);
  }

  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
