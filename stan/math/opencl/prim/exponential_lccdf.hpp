#ifndef STAN_MATH_OPENCL_PRIM_EXPONENTIAL_LCCDF_HPP
#define STAN_MATH_OPENCL_PRIM_EXPONENTIAL_LCCDF_HPP
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
 * Calculates the log exponential cumulative distribution function for
 * the given y and beta.
 *
 * Inverse scale parameter must be greater than 0.
 * y must be greater than or equal to 0.
 *
 * @tparam T_y_cl type of scalar
 * @tparam T_inv_scale_cl type of inverse scale
 * @param y a scalar variable
 * @param beta inverse scale parameter
 * @return The log of the product of densities
 */
template <typename T_y_cl, typename T_inv_scale_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_y_cl, T_inv_scale_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl, T_inv_scale_cl>* = nullptr>
return_type_t<T_y_cl, T_inv_scale_cl> exponential_lccdf(
    const T_y_cl& y, const T_inv_scale_cl& beta) {
  static constexpr const char* function = "exponential_lccdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_inv_scale_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y,
                         "Inverse scale parameter", beta);
  const size_t N = max_size(y, beta);
  if (N == 0) {
    return 0.0;
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

  auto lccdf_expr = colwise_sum(elt_multiply(-beta_val, y_val));

  matrix_cl<double> lccdf_cl;

  results(check_y_nonnegative, check_beta_positive_finite, lccdf_cl)
      = expressions(y_nonnegative_expr, beta_positive_finite_expr, lccdf_expr);

  T_partials_return lccdf = (from_matrix_cl(lccdf_cl)).sum();

  auto ops_partials = make_partials_propagator(y_col, beta_col);

  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = constant(0.0, N, 1) - beta_val;
  }
  if (!is_constant<T_inv_scale_cl>::value) {
    partials<1>(ops_partials) = constant(0.0, N, 1) - y_val;
  }
  return ops_partials.build(lccdf);
}

}  // namespace math
}  // namespace stan
#endif
#endif
