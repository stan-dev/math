#ifndef STAN_MATH_OPENCL_PRIM_UNIFORM_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_UNIFORM_LPDF_HPP
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
 * The log of a uniform density for the given
 * y, lower, and upper bound.
 *
 \f{eqnarray*}{
 y &\sim& \mbox{\sf{U}}(\alpha, \beta) \\
 \log (p (y \, |\, \alpha, \beta)) &=& \log \left( \frac{1}{\beta-\alpha}
 \right) \\
 &=& \log (1) - \log (\beta - \alpha) \\
 &=& -\log (\beta - \alpha) \\
 & & \mathrm{ where } \; y \in [\alpha, \beta], \log(0) \; \mathrm{otherwise}
 \f}
 *
 * @tparam T_y_cl type of scalar
 * @tparam T_low_cl type of lower bound
 * @tparam T_high_cl_cl type of upper bound
 * @param y A scalar variable.
 * @param alpha Lower bound.
 * @param beta Upper bound.
 * @throw std::invalid_argument if the lower bound is greater than
 *    or equal to the lower bound
 */
template <bool propto, typename T_y_cl, typename T_low_cl, typename T_high_cl,
          require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_low_cl,
                                                      T_high_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl, T_low_cl, T_high_cl>* = nullptr>
inline return_type_t<T_y_cl, T_low_cl, T_high_cl> uniform_lpdf(
    const T_y_cl& y, const T_low_cl& alpha, const T_high_cl& beta) {
  static const char* function = "uniform_lpdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_low_cl, T_high_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y,
                         "Lower bound parameter", alpha,
                         "Upper bound parameter", beta);
  const size_t N = max_size(y, alpha, beta);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_low_cl, T_high_cl>::value) {
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
  auto y_not_nan = !isnan(y_val);
  auto check_alpha_finite
      = check_cl(function, "Lower bound parameter", alpha_val, "finite");
  auto alpha_finite = isfinite(alpha_val);
  auto check_beta_finite
      = check_cl(function, "Upper bound parameter", beta_val, "finite");
  auto beta_finite = isfinite(beta_val);

  auto diff = beta_val - alpha_val;

  auto check_diff_positive = check_cl(
      function, "Difference between upper and lower bound", diff, "positive");
  auto diff_positive = diff > 0;

  auto y_out_of_bounds
      = colwise_max(cast<char>(y_val < alpha_val || beta_val < y_val));

  auto logp_expr = colwise_sum(
      static_select<include_summand<propto, T_low_cl, T_high_cl>::value>(
          -log(diff), constant(0.0, N, 1)));

  auto inv_beta_minus_alpha = elt_divide(1.0, diff);

  matrix_cl<char> y_out_of_bounds_cl;
  matrix_cl<double> logp_cl;
  matrix_cl<double> alpha_deriv_cl;
  matrix_cl<double> beta_deriv_cl;

  results(check_y_not_nan, check_alpha_finite, check_beta_finite,
          check_diff_positive, y_out_of_bounds_cl, logp_cl, alpha_deriv_cl,
          beta_deriv_cl)
      = expressions(
          y_not_nan, alpha_finite, beta_finite, diff_positive, y_out_of_bounds,
          logp_expr,
          calc_if<!is_constant<T_low_cl>::value>(inv_beta_minus_alpha),
          calc_if<!is_constant<T_high_cl>::value>(-inv_beta_minus_alpha));

  if (from_matrix_cl(y_out_of_bounds_cl).any()) {
    return LOG_ZERO;
  }

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  operands_and_partials<decltype(y_col), decltype(alpha_col),
                        decltype(beta_col)>
      ops_partials(y_col, alpha_col, beta_col);

  if (!is_constant<T_low_cl>::value) {
    ops_partials.edge2_.partials_ = std::move(alpha_deriv_cl);
  }
  if (!is_constant<T_high_cl>::value) {
    ops_partials.edge3_.partials_ = std::move(beta_deriv_cl);
  }
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
