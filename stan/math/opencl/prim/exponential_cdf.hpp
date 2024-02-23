#ifndef STAN_MATH_OPENCL_PRIM_EXPONENTIAL_CDF_HPP
#define STAN_MATH_OPENCL_PRIM_EXPONENTIAL_CDF_HPP
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
 * Calculates the exponential cumulative distribution function for
 * the given y and beta.
 *
 * Inverse scale parameter must be greater than 0.
 * y must be greater than or equal to 0.
 *
 * @tparam T_y_cl type of scalar
 * @tparam T_inv_scale_cl type of inverse scale
 * @param y a scalar variable
 * @param beta inverse scale parameter
 * @return The product of densities
 */
template <typename T_y_cl, typename T_inv_scale_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_y_cl, T_inv_scale_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl, T_inv_scale_cl>* = nullptr>
return_type_t<T_y_cl, T_inv_scale_cl> exponential_cdf(
    const T_y_cl& y, const T_inv_scale_cl& beta) {
  static constexpr const char* function = "exponential_cdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_inv_scale_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y,
                         "Inverse scale parameter", beta);
  const size_t N = max_size(y, beta);
  if (N == 0) {
    return 1.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& beta_col = as_column_vector_or_scalar(beta);

  const auto& y_val = value_of(y_col);
  const auto& beta_val = value_of(beta_col);

  auto check_y_nonnegative
      = check_cl(function, "Random variable", y_val, "nonnegative");
  auto y_nonnegative_expr = y_val >= 0.0;
  auto check_beta_positive_finite = check_cl(
      function, "Inverse scale parameter", beta_val, "positive finite");
  auto beta_positive_finite_expr = 0.0 < beta_val && isfinite(beta_val);

  auto exp_val = exp(elt_multiply(-beta_val, y_val));
  auto one_m_exp = 1.0 - exp_val;
  auto cdf_expr = colwise_prod(one_m_exp);
  auto rep_deriv1 = elt_divide(exp_val, one_m_exp);

  matrix_cl<double> cdf_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> beta_deriv_cl;

  results(check_y_nonnegative, check_beta_positive_finite, cdf_cl,
          beta_deriv_cl)
      = expressions(
          y_nonnegative_expr, beta_positive_finite_expr, cdf_expr,
          calc_if<!is_constant_all<T_y_cl, T_inv_scale_cl>::value>(rep_deriv1));

  T_partials_return cdf = (from_matrix_cl(cdf_cl)).prod();

  auto rep_deriv2 = beta_deriv_cl * cdf;
  auto y_deriv = elt_multiply(
      static_select<is_constant<T_y_cl>::value>(0, beta_val), rep_deriv2);
  auto beta_deriv = elt_multiply(
      static_select<is_constant<T_inv_scale_cl>::value>(0, y_val), rep_deriv2);

  results(y_deriv_cl, beta_deriv_cl)
      = expressions(calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_inv_scale_cl>::value>(beta_deriv));

  auto ops_partials = make_partials_propagator(y_col, beta_col);

  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  if (!is_constant<T_inv_scale_cl>::value) {
    partials<1>(ops_partials) = std::move(beta_deriv_cl);
  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
#endif
