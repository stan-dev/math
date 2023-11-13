#ifndef STAN_MATH_OPENCL_PRIM_NORMAL_LCCDF_HPP
#define STAN_MATH_OPENCL_PRIM_NORMAL_LCCDF_HPP
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
 * Returns the normal log complementary cumulative distribution function
 * for the given location, and scale. If given containers of matching sizes
 * returns the log sum of probabilities.
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
return_type_t<T_y_cl, T_loc_cl, T_scale_cl> normal_lccdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_scale_cl& sigma) {
  static const char* function = "normal_lccdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_loc_cl, T_scale_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
  const size_t N = max_size(y, mu, sigma);
  if (N == 0) {
    return 0.0;
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
  matrix_cl<double> one_m_erf = select(
      scaled_diff < -37.5 * INV_SQRT_TWO, 2.0,
      select(scaled_diff < -5.0 * INV_SQRT_TWO, 2.0 - erfc(-scaled_diff),
             select(scaled_diff > 8.25 * INV_SQRT_TWO, 0.0,
                    1.0 - erf(scaled_diff))));
  auto lccdf_expr = log(one_m_erf);
  auto mu_deriv = select(scaled_diff > 8.25 * INV_SQRT_TWO, INFTY,
                         SQRT_TWO_OVER_SQRT_PI
                             * elt_divide(exp(-square(scaled_diff)),
                                          elt_multiply(one_m_erf, sigma_val)));
  auto y_deriv = -mu_deriv;
  auto sigma_deriv = elt_multiply(mu_deriv, scaled_diff) * SQRT_TWO;

  matrix_cl<double> lccdf_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;

  results(check_y_not_nan, check_mu_finite, check_sigma_positive)
      = expressions(y_not_nan_expr, mu_finite_expr, sigma_positive_expr);
  results(lccdf_cl, y_deriv_cl, mu_deriv_cl, sigma_deriv_cl)
      = expressions(lccdf_expr, calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_loc_cl>::value>(mu_deriv),
                    calc_if<!is_constant<T_scale_cl>::value>(sigma_deriv));

  T_partials_return lccdf
      = LOG_HALF * lccdf_cl.size() + sum(from_matrix_cl(lccdf_cl));

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
  return ops_partials.build(lccdf);
}

}  // namespace math
}  // namespace stan
#endif
#endif
