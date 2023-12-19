#ifndef STAN_MATH_OPENCL_PRIM_RAYLEIGH_CDF_HPP
#define STAN_MATH_OPENCL_PRIM_RAYLEIGH_CDF_HPP
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
 * Returns the Rayleigh cumulative distribution function for the given
 * location, and scale. If given containers of matching sizes
 * returns the product of probabilities.
 *
 * @tparam T_y_cl type of scalar outcome
 * @tparam T_scale_cl type of scale
 * @param y (Sequence of) scalar(s).
 * @param sigma (Sequence of) scale(s).
 * @return The log of the product of densities.
 */
template <
    typename T_y_cl, typename T_scale_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_scale_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_scale_cl>* = nullptr>
return_type_t<T_y_cl, T_scale_cl> rayleigh_cdf(const T_y_cl& y,
                                               const T_scale_cl& sigma) {
  static const char* function = "rayleigh_cdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_scale_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Scale parameter",
                         sigma);
  const size_t N = max_size(y, sigma);
  if (N == 0) {
    return 1.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& sigma_col = as_column_vector_or_scalar(sigma);

  const auto& y_val = value_of(y_col);
  const auto& sigma_val = value_of(sigma_col);

  auto check_y_nonnegative
      = check_cl(function, "Random variable", y_val, "nonnegative");
  auto y_nonnegative_expr = 0.0 <= y_val;
  auto check_sigma_positive
      = check_cl(function, "Scale parameter", sigma_val, "positive");
  auto sigma_positive_expr = 0 < sigma_val;

  auto inv_sigma = elt_divide(1.0, sigma_val);
  auto inv_sigma_square = square(inv_sigma);
  auto exp_val = exp(elt_multiply(-0.5 * square(y_val), inv_sigma_square));
  auto cdf_expr = colwise_prod(1.0 - exp_val);

  auto y_deriv1 = elt_multiply(elt_multiply(y_val, inv_sigma_square),
                               elt_divide(exp_val, 1.0 - exp_val));
  auto sigma_deriv1 = elt_multiply(elt_multiply(y_val, -inv_sigma), y_deriv1);

  matrix_cl<double> cdf_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;

  results(check_y_nonnegative, check_sigma_positive, cdf_cl, y_deriv_cl,
          sigma_deriv_cl)
      = expressions(y_nonnegative_expr, sigma_positive_expr, cdf_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv1),
                    calc_if<!is_constant<T_scale_cl>::value>(sigma_deriv1));

  T_partials_return cdf = (from_matrix_cl(cdf_cl)).prod();

  auto ops_partials = make_partials_propagator(y_col, sigma_col);

  if (!is_constant_all<T_y_cl, T_scale_cl>::value) {
    results(y_deriv_cl, sigma_deriv_cl) = expressions(
        calc_if<!is_constant<T_y_cl>::value>(y_deriv_cl * cdf),
        calc_if<!is_constant<T_scale_cl>::value>(sigma_deriv_cl * cdf));
    if (!is_constant<T_y_cl>::value) {
      partials<0>(ops_partials) = std::move(y_deriv_cl);
    }
    if (!is_constant<T_scale_cl>::value) {
      partials<1>(ops_partials) = std::move(sigma_deriv_cl);
    }
  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
#endif
