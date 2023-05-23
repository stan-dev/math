#ifndef STAN_MATH_OPENCL_PRIM_STUDENT_T_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_STUDENT_T_LPDF_HPP
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
 * The log of the Student-t density for the given y, nu, mean, and
 * scale parameter.  The scale parameter must be greater
 * than 0.
 *
 * \f{eqnarray*}{
 y &\sim& t_{\nu} (\mu, \sigma^2) \\
 \log (p (y \, |\, \nu, \mu, \sigma) ) &=& \log \left( \frac{\Gamma((\nu + 1)
 /2)}
 {\Gamma(\nu/2)\sqrt{\nu \pi} \sigma} \left( 1 + \frac{1}{\nu} (\frac{y -
 \mu}{\sigma})^2 \right)^{-(\nu + 1)/2} \right) \\
 &=& \log( \Gamma( (\nu+1)/2 )) - \log (\Gamma (\nu/2) - \frac{1}{2} \log(\nu
 \pi) - \log(\sigma)
 -\frac{\nu + 1}{2} \log (1 + \frac{1}{\nu} (\frac{y - \mu}{\sigma})^2)
 \f}
 *
 * @tparam T_y_cl type of scalar
 * @tparam T_dof_cl type of degrees of freedom
 * @tparam T_loc_cl type of location
 * @tparam T_scale_cl type of scale
 *
 * @param y A scalar variable.
 * @param nu Degrees of freedom.
 * @param mu The mean of the Student-t distribution.
 * @param sigma The scale parameter of the Student-t distribution.
 * @return The log of the Student-t density at y.
 * @throw std::domain_error if sigma is not greater than 0.
 * @throw std::domain_error if nu is not greater than 0.
 */
template <bool propto, typename T_y_cl, typename T_dof_cl, typename T_loc_cl,
          typename T_scale_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_y_cl, T_dof_cl, T_loc_cl, T_scale_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl, T_dof_cl, T_loc_cl,
                                        T_scale_cl>* = nullptr>
inline return_type_t<T_y_cl, T_dof_cl, T_loc_cl, T_scale_cl> student_t_lpdf(
    const T_y_cl& y, const T_dof_cl& nu, const T_loc_cl& mu,
    const T_scale_cl& sigma) {
  static const char* function = "student_t_lpdf(OpenCL)";
  using T_partials_return
      = partials_return_t<T_y_cl, T_dof_cl, T_loc_cl, T_scale_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y,
                         "Degrees of freedom parameter", nu,
                         "Location parameter", mu, "Scale parameter", sigma);
  const size_t N = max_size(y, mu, sigma);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_dof_cl, T_loc_cl, T_scale_cl>::value) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& nu_col = as_column_vector_or_scalar(nu);
  const auto& mu_col = as_column_vector_or_scalar(mu);
  const auto& sigma_col = as_column_vector_or_scalar(sigma);

  const auto& y_val = value_of(y_col);
  const auto& nu_val = value_of(nu_col);
  const auto& mu_val = value_of(mu_col);
  const auto& sigma_val = value_of(sigma_col);

  auto check_y_not_nan
      = check_cl(function, "Random variable", y_val, "not NaN");
  auto y_not_nan = !isnan(y_val);
  auto check_nu_positive_finite = check_cl(
      function, "Degrees of freedom parameter", nu_val, "positive finite");
  auto nu_positive_finite = 0 < nu_val && isfinite(nu_val);
  auto check_mu_finite
      = check_cl(function, "Location parameter", mu_val, "finite");
  auto mu_finite = isfinite(mu_val);
  auto check_sigma_positive_finite
      = check_cl(function, "Scale parameter", sigma_val, "positive finite");
  auto sigma_positive_finite = 0 < sigma_val && isfinite(sigma_val);

  auto half_nu = 0.5 * nu_val;
  auto y_scaled = elt_divide(y_val - mu_val, sigma_val);
  auto square_y_scaled = elt_multiply(y_scaled, y_scaled);
  auto square_y_scaled_over_nu = elt_divide(square_y_scaled, nu_val);
  auto log1p_val = log1p(square_y_scaled_over_nu);

  auto logp1 = -elt_multiply((half_nu + 0.5), log1p_val);
  auto logp2 = static_select<include_summand<propto, T_dof_cl>::value>(
      logp1 + lgamma(half_nu + 0.5) - lgamma(half_nu) - 0.5 * log(nu_val),
      logp1);
  auto logp_expr
      = colwise_sum(static_select<include_summand<propto, T_scale_cl>::value>(
          logp2 - log(sigma_val), logp2));

  auto square_sigma = elt_multiply(sigma_val, sigma_val);
  auto deriv_y_mu = elt_divide(
      elt_multiply(nu_val + 1, y_val - mu_val),
      elt_multiply(elt_multiply(1 + square_y_scaled_over_nu, square_sigma),
                   nu_val));
  auto rep_deriv = elt_divide(elt_multiply(nu_val + 1, square_y_scaled_over_nu),
                              1 + square_y_scaled_over_nu)
                   - 1;
  auto nu_deriv = 0.5
                  * (digamma(half_nu + 0.5) - digamma(half_nu) - log1p_val
                     + elt_divide(rep_deriv, nu_val));
  auto sigma_deriv = elt_divide(rep_deriv, sigma_val);

  matrix_cl<double> logp_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> nu_deriv_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;

  results(check_y_not_nan, check_nu_positive_finite, check_mu_finite,
          check_sigma_positive_finite, logp_cl, y_deriv_cl, nu_deriv_cl,
          mu_deriv_cl, sigma_deriv_cl)
      = expressions(y_not_nan, nu_positive_finite, mu_finite,
                    sigma_positive_finite, logp_expr,
                    calc_if<!is_constant<T_y_cl>::value>(-deriv_y_mu),
                    calc_if<!is_constant<T_dof_cl>::value>(nu_deriv),
                    calc_if<!is_constant<T_loc_cl>::value>(deriv_y_mu),
                    calc_if<!is_constant<T_scale_cl>::value>(sigma_deriv));

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  if (include_summand<propto>::value) {
    logp -= LOG_SQRT_PI * N;
  }

  auto ops_partials
      = make_partials_propagator(y_col, nu_col, mu_col, sigma_col);

  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  if (!is_constant<T_dof_cl>::value) {
    partials<1>(ops_partials) = std::move(nu_deriv_cl);
  }
  if (!is_constant<T_loc_cl>::value) {
    partials<2>(ops_partials) = std::move(mu_deriv_cl);
  }
  if (!is_constant<T_scale_cl>::value) {
    partials<3>(ops_partials) = std::move(sigma_deriv_cl);
  }
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
