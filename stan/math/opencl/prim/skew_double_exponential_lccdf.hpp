#ifndef STAN_MATH_OPENCL_PRIM_DOUBLE_SKEW_DOUBLE_EXPONENTIAL_LCCDF_HPP
#define STAN_MATH_OPENCL_PRIM_DOUBLE_SKEW_DOUBLE_EXPONENTIAL_LCCDF_HPP
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
 * Returns the skew double exponential cumulative density function. Given
 * containers of matching sizes, returns the product of probabilities.
 *
 * @tparam T_y_cl type of scalar outcome
 * @tparam T_loc_cl type of location
 * @tparam T_scale_cl type of scale
 * @tparam T_skewness_cl type of inverse scale
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) location(s).
 * @param sigma (Sequence of) scale(s).
 * @param tau (Sequence of) inverse scale(s).
 * @return The log of the product of densities.
 */
template <typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
          typename T_skewness_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_y_cl, T_loc_cl, T_scale_cl, T_skewness_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_scale_cl,
                                        T_skewness_cl>* = nullptr>
return_type_t<T_y_cl, T_loc_cl, T_scale_cl, T_skewness_cl>
skew_double_exponential_lccdf(const T_y_cl& y, const T_loc_cl& mu,
                              const T_scale_cl& sigma,
                              const T_skewness_cl& tau) {
  static const char* function = "skew_double_exponential_lccdf(OpenCL)";
  using T_partials_return
      = partials_return_t<T_y_cl, T_loc_cl, T_scale_cl, T_skewness_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Shape parameter", sigma, "Skewness parameter",
                         tau);
  const size_t N = max_size(y, mu, sigma, tau);
  if (N == 0) {
    return 1.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& mu_col = as_column_vector_or_scalar(mu);
  const auto& sigma_col = as_column_vector_or_scalar(sigma);
  const auto& tau_col = as_column_vector_or_scalar(tau);

  const auto& y_val = value_of(y_col);
  const auto& mu_val = value_of(mu_col);
  const auto& sigma_val = value_of(sigma_col);
  const auto& tau_val = value_of(tau_col);

  auto check_y_not_nan
      = check_cl(function, "Random variable", y_val, "not NaN");
  auto y_not_nan_expr = !isnan(y_val);
  auto check_mu_finite
      = check_cl(function, "Location parameter", mu_val, "finite");
  auto mu_finite_expr = isfinite(mu_val);
  auto check_sigma_positive_finite
      = check_cl(function, "Scale parameter", sigma_val, "positive finite");
  auto sigma_positive_finite_expr = 0.0 < sigma_val && isfinite(sigma_val);
  auto check_tau_bounded = check_cl(function, "Skewness parameter", tau_val,
                                    "in the interval [0, 1]");
  auto tau_bounded_expr = 0.0 < tau_val && tau_val <= 1.0;

  auto inv_sigma = elt_divide(1.0, sigma_val);
  auto y_m_mu = y_val - mu_val;
  auto diff_sign = sign(y_m_mu);
  auto diff_sign_smaller_0 = diff_sign < 0;
  auto abs_diff_y_mu = fabs(y_m_mu);
  auto abs_diff_y_mu_over_sigma = elt_multiply(abs_diff_y_mu, inv_sigma);
  auto expo
      = elt_multiply(diff_sign_smaller_0 + elt_multiply(diff_sign, tau_val),
                     abs_diff_y_mu_over_sigma);
  auto tau_minus_1 = tau_val - 1.0;
  auto inv_exp_2_expo_tau = elt_divide(
      1.0, exp(2.0 * elt_multiply(tau_minus_1, elt_multiply(y_m_mu, inv_sigma)))
               - tau_val);
  auto tau_m1_tau = elt_multiply(tau_minus_1, tau_val);

  auto lccdf_expr = colwise_sum(
      select(y_val <= mu_val, log1m(elt_multiply(tau_val, exp(-2.0 * expo))),
             log1m(tau_val) - 2.0 * expo));

  auto cond = y_val < mu_val;
  auto y_deriv = select(cond,
                        2.0
                            * elt_multiply(elt_multiply(tau_m1_tau, inv_sigma),
                                           inv_exp_2_expo_tau),
                        -2.0 * elt_multiply(inv_sigma, tau_val));
  auto mu_deriv = -y_deriv;
  auto sigma_deriv = select(
      cond,
      -2.0
          * elt_multiply(elt_multiply(tau_m1_tau, y_m_mu),
                         elt_multiply(square(inv_sigma), inv_exp_2_expo_tau)),
      2.0 * elt_multiply(inv_sigma, expo));
  auto tau_deriv = select(
      cond,
      -elt_multiply(
          sigma_val
              + 2.0 * elt_multiply(elt_multiply(tau_val, y_m_mu), diff_sign),
          elt_multiply(inv_sigma, inv_exp_2_expo_tau)),
      elt_divide(1.0, tau_minus_1) - 2.0 * elt_multiply(inv_sigma, y_m_mu));

  matrix_cl<double> lccdf_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;
  matrix_cl<double> tau_deriv_cl;

  results(check_y_not_nan, check_mu_finite, check_sigma_positive_finite,
          check_tau_bounded, lccdf_cl, y_deriv_cl, mu_deriv_cl, sigma_deriv_cl,
          tau_deriv_cl)
      = expressions(y_not_nan_expr, mu_finite_expr, sigma_positive_finite_expr,
                    tau_bounded_expr, lccdf_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_loc_cl>::value>(mu_deriv),
                    calc_if<!is_constant<T_scale_cl>::value>(sigma_deriv),
                    calc_if<!is_constant<T_skewness_cl>::value>(tau_deriv));

  T_partials_return lccdf = (from_matrix_cl(lccdf_cl)).sum();

  operands_and_partials<decltype(y_col), decltype(mu_col), decltype(sigma_col),
                        decltype(tau_col)>
      ops_partials(y_col, mu_col, sigma_col, tau_col);

  if (!is_constant<T_y_cl>::value) {
    ops_partials.edge1_.partials_ = std::move(y_deriv_cl);
  }
  if (!is_constant<T_loc_cl>::value) {
    ops_partials.edge2_.partials_ = std::move(mu_deriv_cl);
  }
  if (!is_constant<T_scale_cl>::value) {
    ops_partials.edge3_.partials_ = std::move(sigma_deriv_cl);
  }
  if (!is_constant<T_skewness_cl>::value) {
    ops_partials.edge4_.partials_ = std::move(tau_deriv_cl);
  }

  return ops_partials.build(lccdf);
}

}  // namespace math
}  // namespace stan
#endif
#endif
