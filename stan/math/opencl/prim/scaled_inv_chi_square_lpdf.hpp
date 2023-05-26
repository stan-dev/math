#ifndef STAN_MATH_OPENCL_PRIM_SCALED_INV_CHI_SQUARE_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_SCALED_INV_CHI_SQUARE_LPDF_HPP
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
 * The log of a scaled inverse chi-squared density for y with the
 * specified degrees of freedom parameter and scale parameter.
 *
 \f{eqnarray*}{
 y &\sim& \mbox{\sf{Inv-}}\chi^2(\nu, s^2) \\
 \log (p (y \, |\, \nu, s)) &=& \log \left( \frac{(\nu / 2)^{\nu / 2}}{\Gamma
 (\nu / 2)} s^\nu y^{- (\nu / 2 + 1)} \exp^{-\nu s^2 / (2y)} \right) \\
 &=& \frac{\nu}{2} \log(\frac{\nu}{2}) - \log (\Gamma (\nu / 2)) + \nu \log(s) -
 (\frac{\nu}{2} + 1) \log(y) - \frac{\nu s^2}{2y} \\ & & \mathrm{ where } \; y >
 0 \f}
 *
 * @tparam T_y_cl type of random variable
 * @tparam T_dof_cl type of degrees of freedom
 * @tparam T_Scale_cl type of scale
 * @param y random variable
 * @param nu degrees of freedom
 * @param s Scale parameter.
 * @throw std::domain_error if nu is not greater than 0
 * @throw std::domain_error if s is not greater than 0.
 * @throw std::domain_error if y is not greater than 0.
 */
template <
    bool propto, typename T_y_cl, typename T_dof_cl, typename T_scale_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_dof_cl,
                                                T_scale_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_dof_cl, T_scale_cl>* = nullptr>
inline return_type_t<T_y_cl, T_dof_cl, T_scale_cl> scaled_inv_chi_square_lpdf(
    const T_y_cl& y, const T_dof_cl& nu, const T_scale_cl& s) {
  static constexpr const char* function = "scaled_inv_chi_square_lpdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_dof_cl, T_scale_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y,
                         "Degrees of freedom parameter", nu, "Scale parameter",
                         s);
  const size_t N = max_size(y, nu, s);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_dof_cl, T_scale_cl>::value) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& nu_col = as_column_vector_or_scalar(nu);
  const auto& s_col = as_column_vector_or_scalar(s);

  const auto& y_val = value_of(y_col);
  const auto& nu_val = value_of(nu_col);
  const auto& s_val = value_of(s_col);

  auto check_y_not_nan
      = check_cl(function, "Random variable", y_val, "not NaN");
  auto y_not_nan = !isnan(y_val);
  auto check_nu_positive_finite = check_cl(
      function, "Degrees of freedom parameter", nu_val, "positive finite");
  auto nu_positive_finite = isfinite(nu_val) && 0 < nu_val;
  auto check_s_positive_finite
      = check_cl(function, "Scale parameter", s_val, "positive finite");
  auto s_positive_finite = isfinite(s_val) && 0 < s_val;

  auto any_y_nonpositive = colwise_max(cast<char>(y_val <= 0.0));
  auto half_nu = 0.5 * nu_val;
  auto log_y = log(y_val);
  auto inv_y = elt_divide(1.0, y_val);
  auto log_s = log(s_val);
  auto log_half_nu = log(half_nu);
  auto s_square = elt_multiply(s_val, s_val);

  auto logp1 = static_select<include_summand<propto, T_dof_cl>::value>(
      elt_multiply(half_nu, log_half_nu) - lgamma(half_nu), constant(0, N, 1));
  auto logp2
      = static_select<include_summand<propto, T_dof_cl, T_scale_cl>::value>(
          logp1 + elt_multiply(nu_val, log_s), logp1);
  auto logp3 = static_select<include_summand<propto, T_dof_cl, T_y_cl>::value>(
      logp2 - elt_multiply(half_nu + 1.0, log_y), logp2);
  auto logp_expr = colwise_sum(
      static_select<
          include_summand<propto, T_dof_cl, T_y_cl, T_scale_cl>::value>(
          logp3 - elt_multiply(elt_multiply(half_nu, s_square), inv_y), logp3));

  auto y_deriv = -elt_multiply(half_nu + 1.0, inv_y)
                 + elt_multiply(elt_multiply(half_nu, s_square),
                                elt_multiply(inv_y, inv_y));
  auto nu_deriv = 0.5
                      * (log_half_nu - digamma(half_nu) - log_y
                         - elt_multiply(s_square, inv_y))
                  + log_s + 0.5;
  auto s_deriv = elt_divide(nu_val, s_val)
                 - elt_multiply(elt_multiply(nu_val, inv_y), s_val);

  matrix_cl<char> any_y_nonpositive_cl;
  matrix_cl<double> logp_cl;
  matrix_cl<double> nu_deriv_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> s_deriv_cl;

  results(check_y_not_nan, check_nu_positive_finite, check_s_positive_finite,
          any_y_nonpositive_cl, logp_cl, y_deriv_cl, nu_deriv_cl, s_deriv_cl)
      = expressions(y_not_nan, nu_positive_finite, s_positive_finite,
                    any_y_nonpositive, logp_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_dof_cl>::value>(nu_deriv),
                    calc_if<!is_constant<T_scale_cl>::value>(s_deriv));

  if (from_matrix_cl(any_y_nonpositive_cl).any()) {
    return LOG_ZERO;
  }

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  auto ops_partials = make_partials_propagator(y_col, nu_col, s_col);

  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  if (!is_constant<T_dof_cl>::value) {
    partials<1>(ops_partials) = std::move(nu_deriv_cl);
  }
  if (!is_constant<T_scale_cl>::value) {
    partials<2>(ops_partials) = std::move(s_deriv_cl);
  }
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
