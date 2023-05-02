#ifndef STAN_MATH_OPENCL_PRIM_LOGNORMAL_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_LOGNORMAL_LPDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * The log of the lognormal density for the specified scalar(s) given the
 * specified sample stan::math::size(s). y, mu, or sigma can each either be
 * scalar or a vector on OpenCL device. Any vector inputs must be the same
 * length.
 *
 * <p> The result log probability is defined to be the sum of
 * the log probabilities for each observation/mu/sigma triple.
 *
 * @tparam T_y_cl type of scalar outcome
 * @tparam T_loc_cl type of prior scale for successes
 * @tparam T_scale_cl type of prior scale for failures
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) prior sample stan::math::size(s).
 * @param sigma (Sequence of) prior sample stan::math::size(s).
 * @return The log of the product of densities.
 */
template <
    bool propto, typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_loc_cl,
                                                T_scale_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_scale_cl>* = nullptr>
return_type_t<T_y_cl, T_loc_cl, T_scale_cl> lognormal_lpdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_scale_cl& sigma) {
  using std::isfinite;
  static const char* function = "lognormal_lpdf(OpenCL)";
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

  auto ops_partials = make_partials_propagator(y_col, mu_col, sigma_col);

  auto check_y_nonnegative
      = check_cl(function, "Random variable", y_val, "nonnegative");
  auto y_nonnegative = 0 <= y_val;
  auto check_mu_finite
      = check_cl(function, "Location parameter", mu_val, "finite");
  auto mu_finite = isfinite(mu_val);
  auto check_sigma_pos_finite
      = check_cl(function, "Scale parameter", sigma_val, "positive finite");
  auto sigma_pos_finite = sigma_val > 0 && isfinite(sigma_val);

  auto any_y_zero = colwise_max(cast<char>(y_val == 0.0));
  auto inv_sigma = elt_divide(1.0, sigma_val);
  auto inv_sigma_sq = elt_multiply(inv_sigma, inv_sigma);
  auto log_y = log(y_val);
  auto logy_m_mu = log_y - mu_val;
  auto logy_m_mu_div_sigma = elt_multiply(logy_m_mu, inv_sigma_sq);

  auto logp1
      = -0.5 * elt_multiply(elt_multiply(logy_m_mu, logy_m_mu), inv_sigma_sq);
  auto logp2 = static_select<include_summand<propto, T_scale_cl>::value>(
      logp1 - log(sigma_val), logp1);
  auto logp_expr
      = colwise_sum(static_select<include_summand<propto, T_y_cl>::value>(
          logp2 - log_y, logp2));

  auto y_deriv_expr = elt_divide(-(1.0 + logy_m_mu_div_sigma), y_val);
  auto sigma_deriv_expr = elt_multiply(
      elt_multiply(logy_m_mu_div_sigma, logy_m_mu) - 1.0, inv_sigma);

  matrix_cl<char> any_y_zero_cl;
  matrix_cl<double> logp_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;

  results(check_y_nonnegative, check_mu_finite, check_sigma_pos_finite,
          any_y_zero_cl, logp_cl, y_deriv_cl, mu_deriv_cl, sigma_deriv_cl)
      = expressions(y_nonnegative, mu_finite, sigma_pos_finite, any_y_zero,
                    logp_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv_expr),
                    calc_if<!is_constant<T_loc_cl>::value>(logy_m_mu_div_sigma),
                    calc_if<!is_constant<T_scale_cl>::value>(sigma_deriv_expr));

  if (from_matrix_cl(any_y_zero_cl).any()) {
    return LOG_ZERO;
  }

  T_partials_return logp
      = sum(from_matrix_cl(logp_cl)) + N * NEG_LOG_SQRT_TWO_PI;

  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  if (!is_constant<T_loc_cl>::value) {
    partials<1>(ops_partials) = std::move(mu_deriv_cl);
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
