#ifndef STAN_MATH_OPENCL_PRIM_GUMBEL_LCCDF_HPP
#define STAN_MATH_OPENCL_PRIM_GUMBEL_LCCDF_HPP
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
 * Returns the Gumbel log complementary cumulative distribution function for the
 * given location, and scale. If given containers of matching sizes returns the
 * product of probabilities.
 *
 * @tparam T_y_cl type of scalar outcome
 * @tparam T_loc_cl type of location
 * @tparam T_scale_cl type of scale
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) location(s).
 * @param beta (Sequence of) scale(s).
 * @return The sum of log complementary cumulative probabilities.
 */
template <
    typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_loc_cl,
                                                T_scale_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_scale_cl>* = nullptr>
return_type_t<T_y_cl, T_loc_cl, T_scale_cl> gumbel_lccdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_scale_cl& beta) {
  static const char* function = "gumbel_lccdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_loc_cl, T_scale_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", beta);
  const size_t N = max_size(y, mu, beta);
  if (N == 0) {
    return 1.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& mu_col = as_column_vector_or_scalar(mu);
  const auto& beta_col = as_column_vector_or_scalar(beta);

  const auto& y_val = value_of(y_col);
  const auto& mu_val = value_of(mu_col);
  const auto& beta_val = value_of(beta_col);

  auto check_y_not_nan
      = check_cl(function, "Random variable", y_val, "not NaN");
  auto y_not_nan_expr = !isnan(y_val);
  auto check_mu_finite
      = check_cl(function, "Location parameter", mu_val, "finite");
  auto mu_finite_expr = isfinite(mu_val);
  auto check_beta_positive
      = check_cl(function, "Scale parameter", beta_val, "positive");
  auto beta_positive_expr = 0.0 < beta_val;

  auto scaled_diff = elt_divide(y_val - mu_val, beta_val);
  auto exp_m_scaled_diff = exp(-scaled_diff);
  auto ccdf_n = 1.0 - exp(-exp_m_scaled_diff);
  auto lccdf_expr = colwise_sum(log(ccdf_n));
  auto mu_deriv = elt_divide(exp(-scaled_diff - exp_m_scaled_diff),
                             elt_multiply(beta_val, ccdf_n));
  auto y_deriv = -mu_deriv;
  auto beta_deriv = elt_multiply(mu_deriv, scaled_diff);

  matrix_cl<double> lccdf_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> beta_deriv_cl;

  results(check_y_not_nan, check_mu_finite, check_beta_positive, lccdf_cl,
          y_deriv_cl, mu_deriv_cl, beta_deriv_cl)
      = expressions(y_not_nan_expr, mu_finite_expr, beta_positive_expr,
                    lccdf_expr, calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_loc_cl>::value>(mu_deriv),
                    calc_if<!is_constant<T_scale_cl>::value>(beta_deriv));

  T_partials_return lccdf = sum(from_matrix_cl(lccdf_cl));

  auto ops_partials = make_partials_propagator(y_col, mu_col, beta_col);

  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  if (!is_constant<T_loc_cl>::value) {
    partials<1>(ops_partials) = std::move(mu_deriv_cl);
  }
  if (!is_constant<T_scale_cl>::value) {
    partials<2>(ops_partials) = std::move(beta_deriv_cl);
  }
  return ops_partials.build(lccdf);
}

}  // namespace math
}  // namespace stan
#endif
#endif
