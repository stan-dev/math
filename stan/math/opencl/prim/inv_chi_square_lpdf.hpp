#ifndef STAN_MATH_OPENCL_PRIM_INV_CHI_SQUARE_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_INV_CHI_SQUARE_LPDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * The log of an inverse chi-squared density for y with the specified
 * degrees of freedom parameter.
 * The degrees of freedom parameter must be greater than 0.
 * y must be greater than 0.
 *
 \f{eqnarray*}{
 y &\sim& \mbox{\sf{Inv-}}\chi^2_\nu \\
 \log (p (y \, |\, \nu)) &=& \log \left( \frac{2^{-\nu / 2}}{\Gamma (\nu / 2)}
 y^{- (\nu / 2 + 1)} \exp^{-1 / (2y)} \right) \\
 &=& - \frac{\nu}{2} \log(2) - \log (\Gamma (\nu / 2)) - (\frac{\nu}{2} + 1)
 \log(y) - \frac{1}{2y} \\ & & \mathrm{ where } \; y > 0 \f}
 *
 * @tparam T_y type of scalar
 * @tparam T_dof_cl type of degrees of freedom
 * @param y A scalar variable.
 * @param nu Degrees of freedom.
 * @throw std::domain_error if nu is not greater than or equal to 0
 * @throw std::domain_error if y is not greater than or equal to 0.
 */
template <
    bool propto, typename T_y_cl, typename T_dof_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_dof_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_dof_cl>* = nullptr>
return_type_t<T_y_cl, T_dof_cl> inv_chi_square_lpdf(const T_y_cl& y,
                                                    const T_dof_cl& nu) {
  using std::isfinite;
  using std::isnan;
  static const char* function = "inv_chi_square_lpdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_dof_cl>;

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

  auto ops_partials = make_partials_propagator(y_col, nu_col);

  auto check_nu_pos_finite = check_cl(function, "Degrees of freedom parameter",
                                      nu_val, "positive finite");
  auto nu_pos_finite = nu_val > 0 && isfinite(nu_val);
  auto check_y_not_nan
      = check_cl(function, "Random variable", y_val, "not NaN");
  auto y_not_nan = !isnan(y_val);

  auto any_y_nonpositive = colwise_max(cast<char>(y_val <= 0));
  auto log_y = log(y_val);
  auto half_nu = nu_val * 0.5;
  auto two_over_y = elt_divide(0.5, y_val);

  auto logp1 = -elt_multiply(half_nu + 1.0, log_y);
  auto logp2 = static_select<include_summand<propto, T_dof_cl>::value>(
      logp1 - nu_val * HALF_LOG_TWO - lgamma(half_nu), logp1);
  auto logp_expr
      = colwise_sum(static_select<include_summand<propto, T_y_cl>::value>(
          logp2 - two_over_y, logp2));

  auto y_deriv = elt_divide(two_over_y - half_nu - 1.0, y_val);
  auto nu_deriv = -HALF_LOG_TWO - (digamma(half_nu) + log_y) * 0.5;

  matrix_cl<char> any_y_nonpositive_cl;
  matrix_cl<double> logp_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> nu_deriv_cl;

  results(check_nu_pos_finite, check_y_not_nan, any_y_nonpositive_cl, logp_cl,
          y_deriv_cl, nu_deriv_cl)
      = expressions(nu_pos_finite, y_not_nan, any_y_nonpositive, logp_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_dof_cl>::value>(nu_deriv));

  if (from_matrix_cl(any_y_nonpositive_cl).any()) {
    return LOG_ZERO;
  }

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  if (!is_constant<T_dof_cl>::value) {
    partials<1>(ops_partials) = std::move(nu_deriv_cl);
  }

  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
