#ifndef STAN_MATH_OPENCL_PRIM_NORMAL_CDF_HPP
#define STAN_MATH_OPENCL_PRIM_NORMAL_CDF_HPP
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
 * Returns the normal cumulative distribution function for the given
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
return_type_t<T_y_cl, T_loc_cl, T_scale_cl> normal_cdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_scale_cl& sigma) {
  static const char* function = "normal_cdf(OpenCL)";
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

  auto check_y_not_nan
      = check_cl(function, "Random variable", y_val, "not NaN");
  auto y_not_nan_expr = !isnan(y_val);
  auto check_mu_finite
      = check_cl(function, "Location parameter", mu_val, "finite");
  auto mu_finite_expr = isfinite(mu_val);
  auto check_sigma_positive
      = check_cl(function, "Scale parameter", sigma_val, "positive");
  auto sigma_positive_expr = 0 < sigma_val;

  auto scaled_diff = elt_divide(y_val - mu_val, sigma_val * SQRT_TWO);
  auto cdf_n = select(
      scaled_diff < -37.5 * INV_SQRT_TWO, 0.0,
      select(scaled_diff < -5.0 * INV_SQRT_TWO, 0.5 * erfc(-scaled_diff),
             select(scaled_diff > 8.25 * INV_SQRT_TWO, 1.0,
                    0.5 * (1.0 + erf(scaled_diff)))));
  auto cdf_expr = colwise_prod(cdf_n);
  auto mu_deriv_tmp = select(scaled_diff < -37.5 * INV_SQRT_TWO, 0.0,
                             INV_SQRT_TWO_PI
                                 * elt_divide(exp(-square(scaled_diff)),
                                              elt_multiply(cdf_n, sigma_val)));
  auto sigma_deriv_tmp = elt_multiply(mu_deriv_tmp, scaled_diff);

  matrix_cl<double> cdf_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;

  results(check_y_not_nan, check_mu_finite, check_sigma_positive, cdf_cl,
          mu_deriv_cl, sigma_deriv_cl)
      = expressions(
          y_not_nan_expr, mu_finite_expr, sigma_positive_expr, cdf_expr,
          calc_if<!is_constant_all<T_y_cl, T_loc_cl>::value>(mu_deriv_tmp),
          calc_if<!is_constant<T_scale_cl>::value>(sigma_deriv_tmp));

  T_partials_return cdf = (from_matrix_cl(cdf_cl)).prod();

  auto y_deriv = elt_multiply(mu_deriv_cl, cdf);
  auto mu_deriv = -y_deriv;
  auto sigma_deriv = elt_multiply(sigma_deriv_cl, -SQRT_TWO * cdf);

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
