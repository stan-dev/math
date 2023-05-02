#ifndef STAN_MATH_OPENCL_PRIM_STD_NORMAL_CDF_HPP
#define STAN_MATH_OPENCL_PRIM_STD_NORMAL_CDF_HPP
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
 * Returns the standard normal cumulative distribution function.
 *
 * @tparam T_y_cl type of scalar outcome
 * @param y (Sequence of) scalar(s).
 * @return The log of the product of densities.
 */
template <typename T_y_cl,
          require_all_prim_or_rev_kernel_expression_t<T_y_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl>* = nullptr>
return_type_t<T_y_cl> std_normal_cdf(const T_y_cl& y) {
  static const char* function = "std_normal_cdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl>;
  using std::isfinite;
  using std::isnan;

  const size_t N = math::size(y);
  if (N == 0) {
    return 1.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& y_val = value_of(y_col);

  auto check_y_not_nan
      = check_cl(function, "Random variable", y_val, "not NaN");
  auto y_not_nan_expr = !isnan(y_val);

  auto scaled_y = y_val * INV_SQRT_TWO;
  auto cdf_n
      = select(y_val < -37.5, 0.0,
               select(y_val < -5.0, 0.5 * erfc(-scaled_y),
                      select(y_val > 8.25, 1.0, 0.5 * (1.0 + erf(scaled_y)))));
  auto cdf_expr = colwise_prod(cdf_n);
  auto y_deriv1
      = select(y_val < -37.5, 0.0,
               INV_SQRT_TWO_PI * elt_divide(exp(-square(scaled_y)), cdf_n));

  matrix_cl<double> cdf_cl;
  matrix_cl<double> y_deriv_cl;

  results(check_y_not_nan, cdf_cl, y_deriv_cl) = expressions(
      y_not_nan_expr, cdf_expr, calc_if<!is_constant<T_y_cl>::value>(y_deriv1));

  T_partials_return cdf = (from_matrix_cl(cdf_cl)).prod();

  auto ops_partials = make_partials_propagator(y_col);

  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = y_deriv_cl * cdf;
  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
#endif
