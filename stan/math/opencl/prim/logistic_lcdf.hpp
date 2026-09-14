#ifndef STAN_MATH_OPENCL_PRIM_LOGISTIC_LCDF_HPP
#define STAN_MATH_OPENCL_PRIM_LOGISTIC_LCDF_HPP
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
 * Returns the logistic cumulative distribution function for the given
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
inline return_type_t<T_y_cl, T_loc_cl, T_scale_cl> logistic_lcdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_scale_cl& sigma) {
  static constexpr const char* function = "logistic_lcdf(OpenCL)";
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
  auto check_sigma_positive_finite
      = check_cl(function, "Scale parameter", sigma_val, "positive finite");
  auto sigma_positive_finite_expr = 0 < sigma_val && isfinite(sigma_val);

  auto any_y_neg_inf = colwise_max(cast<char>(y_val == NEGATIVE_INFTY));
  auto cond = y_val == INFTY;
  auto inv_sigma = elt_divide(1.0, sigma_val);
  auto scaled_diff = elt_multiply(y_val - mu_val, inv_sigma);
  auto P_expr = colwise_sum(log_inv_logit(scaled_diff));

  // y == INFTY contributes log(1) = 0 to P and zero to every partial; without
  // the select the scale partial would be 0 * INFTY = NaN, which prim (where
  // the element is skipped outright) never produces.
  auto deriv = elt_multiply(inv_logit(-scaled_diff), inv_sigma);
  auto y_deriv = select(cond, 0.0, deriv);
  auto mu_deriv = -y_deriv;
  auto sigma_deriv = select(cond, 0.0, elt_multiply(-deriv, scaled_diff));

  matrix_cl<char> any_y_neg_inf_cl;
  matrix_cl<double> P_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;

  results(check_y_not_nan, check_mu_finite, check_sigma_positive_finite,
          any_y_neg_inf_cl, P_cl, y_deriv_cl, mu_deriv_cl, sigma_deriv_cl)
      = expressions(y_not_nan_expr, mu_finite_expr, sigma_positive_finite_expr,
                    any_y_neg_inf, P_expr,
                    calc_if<is_autodiff_v<T_y_cl>>(y_deriv),
                    calc_if<is_autodiff_v<T_loc_cl>>(mu_deriv),
                    calc_if<is_autodiff_v<T_scale_cl>>(sigma_deriv));

  if (from_matrix_cl(any_y_neg_inf_cl).maxCoeff()) {
    return NEGATIVE_INFTY;
  }

  T_partials_return P = (from_matrix_cl(P_cl)).sum();

  auto ops_partials = make_partials_propagator(y_col, mu_col, sigma_col);

  if constexpr (is_autodiff_v<T_y_cl>) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  if constexpr (is_autodiff_v<T_loc_cl>) {
    partials<1>(ops_partials) = std::move(mu_deriv_cl);
  }
  if constexpr (is_autodiff_v<T_scale_cl>) {
    partials<2>(ops_partials) = std::move(sigma_deriv_cl);
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
#endif
