#ifndef STAN_MATH_OPENCL_PRIM_STD_NORMAL_LCDF_HPP
#define STAN_MATH_OPENCL_PRIM_STD_NORMAL_LCDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/kernels/device_functions/std_normal_lcdf.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {
namespace internal {
constexpr char std_normal_lcdf_opencl_func[] = "std_normal_lcdf(OpenCL)";
}  // namespace internal
/** \ingroup opencl
 * Returns the log standard normal complementary cumulative distribution
 * function.
 *
 * @tparam T_y_cl type of scalar outcome
 * @param y (Sequence of) scalar(s).
 * @return The log of the product of densities.
 */
template <const char* func = internal::std_normal_lcdf_opencl_func,
          typename T_y_cl,
          require_all_prim_or_rev_kernel_expression_t<T_y_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl>* = nullptr>
inline return_type_t<T_y_cl> std_normal_lcdf(const T_y_cl& y) {
  static constexpr const char* function = func;
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
  auto lcdf_expr = colwise_sum(std_normal_lcdf_scaled_impl(scaled_y));
  auto dnlcdf = std_normal_lcdf_dscaled_impl(scaled_y);
  auto y_deriv = dnlcdf * INV_SQRT_TWO;

  matrix_cl<double> lcdf_cl;
  matrix_cl<double> y_deriv_cl;

  results(check_y_not_nan, lcdf_cl, y_deriv_cl) = expressions(
      y_not_nan_expr, lcdf_expr, calc_if<is_autodiff_v<T_y_cl>>(y_deriv));

  double lcdf = from_matrix_cl(lcdf_cl).sum();

  auto ops_partials = make_partials_propagator(y_col);

  if constexpr (is_autodiff_v<T_y_cl>) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  return ops_partials.build(lcdf);
}

}  // namespace math
}  // namespace stan
#endif
#endif
