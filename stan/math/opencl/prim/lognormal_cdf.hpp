#ifndef STAN_MATH_OPENCL_PRIM_LOGNORMAL_CDF_HPP
#define STAN_MATH_OPENCL_PRIM_LOGNORMAL_CDF_HPP
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
 * Returns the loghormal cumulative distribution function for the given
 * location, and scale. If given containers of matching sizes
 * returns the product of probabilities.
 *
 * @tparam T_y_cl type of scalar outcome
 * @tparam T_loc_cl type of location
 * @tparam T_scale_cl type of scale
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) location(s).
 * @param sigma (Sequence of) scale(s).
 * @return The log of the product of densities.
 */
template <
    typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_loc_cl,
                                                T_scale_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_scale_cl>* = nullptr>
return_type_t<T_y_cl, T_loc_cl, T_scale_cl> lognormal_cdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_scale_cl& sigma) {
  static constexpr const char* function = "lognormal_cdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_loc_cl, T_scale_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
  const size_t N = max_size(y, mu, sigma);
  if (N == 0) {
    return 1.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& mu_col = as_column_vector_or_scalar(mu);
  const auto& sigma_col = as_column_vector_or_scalar(sigma);

  const auto& y_val = value_of(y_col);
  const auto& mu_val = value_of(mu_col);
  const auto& sigma_val = value_of(sigma_col);

  auto check_y_nonnegative
      = check_cl(function, "Random variable", y_val, "nonnegative");
  auto y_nonnegative_expr = 0.0 <= y_val;
  auto check_mu_finite
      = check_cl(function, "Location parameter", mu_val, "finite");
  auto mu_finite_expr = isfinite(mu_val);
  auto check_sigma_positive_finite
      = check_cl(function, "Scale parameter", sigma_val, "positive finite");
  auto sigma_positive_finite_expr = 0 < sigma_val && isfinite(sigma_val);

  auto any_y_zero = colwise_max(cast<char>(y_val == 0.0));
  auto log_y = log(y_val);
  auto scaled_diff = elt_divide(log_y - mu_val, sigma_val * SQRT_TWO);
  auto erfc_m_diff = erfc(-scaled_diff);
  auto cdf_n = 0.5 * erfc_m_diff;
  auto cdf_expr = colwise_prod(cdf_n);
  auto mu_deriv_tmp
      = -INV_SQRT_TWO_PI
        * elt_divide(exp(-square(scaled_diff)), elt_multiply(sigma_val, cdf_n));
  auto y_deriv_tmp = elt_divide(-mu_deriv_tmp, y_val);
  auto sigma_deriv_tmp = elt_multiply(mu_deriv_tmp, scaled_diff * SQRT_TWO);

  matrix_cl<char> any_y_zero_cl;
  matrix_cl<double> cdf_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;

  results(check_y_nonnegative, check_mu_finite, check_sigma_positive_finite,
          any_y_zero_cl, cdf_cl, y_deriv_cl, mu_deriv_cl, sigma_deriv_cl)
      = expressions(y_nonnegative_expr, mu_finite_expr,
                    sigma_positive_finite_expr, any_y_zero, cdf_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv_tmp),
                    calc_if<!is_constant<T_loc_cl>::value>(mu_deriv_tmp),
                    calc_if<!is_constant<T_scale_cl>::value>(sigma_deriv_tmp));

  if (from_matrix_cl(any_y_zero_cl).maxCoeff()) {
    return 0.0;
  }

  T_partials_return cdf = (from_matrix_cl(cdf_cl)).prod();

  auto mu_deriv = mu_deriv_cl * cdf;
  auto y_deriv = y_deriv_cl * cdf;
  auto sigma_deriv = sigma_deriv_cl * cdf;

  results(mu_deriv_cl, y_deriv_cl, sigma_deriv_cl)
      = expressions(calc_if<!is_constant<T_loc_cl>::value>(mu_deriv),
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_scale_cl>::value>(sigma_deriv));

  auto ops_partials = make_partials_propagator(y_col, mu_col, sigma_col);

  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  if (!is_constant<T_loc_cl>::value) {
    partials<1>(ops_partials) = std::move(mu_deriv_cl);
  }
  if (!is_constant<T_scale_cl>::value) {
    partials<2>(ops_partials) = std::move(sigma_deriv_cl);
  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
#endif
