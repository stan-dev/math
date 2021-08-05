#ifndef STAN_MATH_OPENCL_PRIM_CHI_SQUARE_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_CHI_SQUARE_LPDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * The log of a chi-squared density for y with the specified
 * degrees of freedom parameter.
 * The degrees of freedom parameter must be greater than 0.
 * y must be greater than or equal to 0.
 *
 \f{eqnarray*}{
 y &\sim& \chi^2_\nu \\
 \log (p (y \, |\, \nu)) &=& \log \left( \frac{2^{-\nu / 2}}{\Gamma (\nu / 2)}
 y^{\nu / 2 - 1} \exp^{- y / 2} \right) \\
 &=& - \frac{\nu}{2} \log(2) - \log (\Gamma (\nu / 2)) + (\frac{\nu}{2} - 1)
 \log(y) - \frac{y}{2} \\ & & \mathrm{ where } \; y \ge 0 \f}
 *
 * @tparam T_y_cl type of dependent variable
 * @tparam T_dof_cl type of degrees of freedom
 * @param y A dependent variable.
 * @param nu Degrees of freedom.
 * @throw std::domain_error if nu is not greater than or equal to 0
 * @throw std::domain_error if y is not greater than or equal to 0.
 */
template <
    bool propto, typename T_y_cl, typename T_dof_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_dof_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_dof_cl>* = nullptr>
return_type_t<T_y_cl, T_dof_cl> chi_square_lpdf(const T_y_cl& y,
                                                const T_dof_cl& nu) {
  static const char* function = "chi_square_lpdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_dof_cl>;
  using std::isfinite;

  check_consistent_sizes(function, "Random variable", y,
                         "Degrees of freedom parameter", nu);
  const size_t N = max_size(y, nu);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_dof_cl>::value) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& nu_col = as_column_vector_or_scalar(nu);

  const auto& y_val = value_of(y_col);
  const auto& nu_val = value_of(nu_col);

  auto check_y_nonnegative
      = check_cl(function, "Random variable", y_val, "nonnegative");
  auto y_nonnegative_expr = 0 <= y_val;
  auto check_nu_positive_finite = check_cl(
      function, "Degrees of freedom parameter", nu_val, "positive finite");
  auto nu_positive_finite_expr = 0 < nu_val && isfinite(nu_val);

  auto log_y_expr = log(y_val);
  auto half_nu_expr = 0.5 * nu_val;

  auto logp1_expr = elt_multiply(half_nu_expr - 1, log_y_expr);
  auto logp2_expr = static_select<include_summand<propto, T_dof_cl>::value>(
      logp1_expr - nu_val * HALF_LOG_TWO - lgamma(half_nu_expr), logp1_expr);
  auto logp_expr
      = colwise_sum(static_select<include_summand<propto, T_y_cl>::value>(
          logp2_expr - 0.5 * y_val, logp2_expr));

  auto y_deriv_expr = elt_divide(half_nu_expr - 1, y_val) - 0.5;
  auto nu_deriv_expr
      = (log_y_expr - digamma(half_nu_expr)) * 0.5 - HALF_LOG_TWO;

  matrix_cl<double> logp_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> nu_deriv_cl;

  results(check_y_nonnegative, check_nu_positive_finite, logp_cl, y_deriv_cl,
          nu_deriv_cl)
      = expressions(y_nonnegative_expr, nu_positive_finite_expr, logp_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv_expr),
                    calc_if<!is_constant<T_dof_cl>::value>(nu_deriv_expr));

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  operands_and_partials<decltype(y_col), decltype(nu_col)> ops_partials(y_col,
                                                                        nu_col);

  if (!is_constant<T_y_cl>::value) {
    ops_partials.edge1_.partials_ = std::move(y_deriv_cl);
  }
  if (!is_constant<T_dof_cl>::value) {
    ops_partials.edge2_.partials_ = std::move(nu_deriv_cl);
  }
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan

#endif
#endif
