#ifndef STAN_MATH_OPENCL_PRIM_EXPONENTIAL_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_EXPONENTIAL_LPDF_HPP
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
 * The log of an exponential density for y with the specified
 * inverse scale parameter.
 * Inverse scale parameter must be greater than 0.
 * y must be greater than or equal to 0.
 *
 \f{eqnarray*}{
 y
 &\sim&
 \mbox{\sf{Expon}}(\beta) \\
 \log (p (y \, |\, \beta) )
 &=&
 \log \left( \beta \exp^{-\beta y} \right) \\
 &=&
 \log (\beta) - \beta y \\
 & &
 \mathrm{where} \; y > 0
 \f}
 *
 * @tparam T_y_cl type of scalar
 * @tparam T_inv_scale_cl type of inverse scale
 * @param y A scalar variable.
 * @param beta Inverse scale parameter.
 * @throw std::domain_error if beta is not greater than 0.
 * @throw std::domain_error if y is not greater than or equal to 0.
 */
template <bool propto, typename T_y_cl, typename T_inv_scale_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_y_cl, T_inv_scale_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl, T_inv_scale_cl>* = nullptr>
return_type_t<T_y_cl, T_inv_scale_cl> exponential_lpdf(
    const T_y_cl& y, const T_inv_scale_cl& beta) {
  using std::isfinite;
  static constexpr const char* function = "exponential_lpdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_inv_scale_cl>;

  check_consistent_sizes(function, "Random variable", y,
                         "Inverse scale parameter", beta);
  const size_t N = max_size(y, beta);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_inv_scale_cl>::value) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& beta_col = as_column_vector_or_scalar(beta);

  const auto& y_val = value_of(y_col);
  const auto& beta_val = value_of(beta_col);

  auto ops_partials = make_partials_propagator(y_col, beta_col);

  auto check_y_nonnegative
      = check_cl(function, "Random variable", y_val, "nonnegative");
  auto y_nonnegative_expr = y_val >= 0;
  auto check_beta_pos_finite = check_cl(function, "Inverse scale parameter",
                                        beta_val, "positive finite");
  auto beta_pos_finite_expr = beta_val > 0 && isfinite(beta_val);

  auto logp1_expr
      = static_select<include_summand<propto, T_inv_scale_cl>::value>(
          log(beta_val), 0);
  auto logp_expr = colwise_sum(
      static_select<include_summand<propto, T_y_cl, T_inv_scale_cl>::value>(
          logp1_expr - elt_multiply(beta_val, y_val), logp1_expr));

  auto y_deriv_expr = elt_multiply(beta_val, constant(-1.0, N, 1));
  auto beta_deriv_expr = elt_divide(1.0, beta_val) - y_val;

  matrix_cl<double> logp_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> beta_deriv_cl;

  results(check_y_nonnegative, check_beta_pos_finite, logp_cl, y_deriv_cl,
          beta_deriv_cl)
      = expressions(
          y_nonnegative_expr, beta_pos_finite_expr, logp_expr,
          calc_if<!is_constant<T_y_cl>::value>(y_deriv_expr),
          calc_if<!is_constant<T_inv_scale_cl>::value>(beta_deriv_expr));

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  if (!is_constant<T_inv_scale_cl>::value) {
    partials<1>(ops_partials) = std::move(beta_deriv_cl);
  }

  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
