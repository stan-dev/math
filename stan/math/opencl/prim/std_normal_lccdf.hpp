#ifndef STAN_MATH_OPENCL_PRIM_STD_NORMAL_LCCDF_HPP
#define STAN_MATH_OPENCL_PRIM_STD_NORMAL_LCCDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns the log standard normal complementary cumulative distribution
 * function.
 *
 * @tparam T_y_cl type of scalar outcome
 * @param y (Sequence of) scalar(s).
 * @return The log of the product of densities.
 */
template <typename T_y_cl,
          require_all_prim_or_rev_kernel_expression_t<T_y_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl>* = nullptr>
return_type_t<T_y_cl> std_normal_lccdf(const T_y_cl& y) {
  static const char* function = "std_normal_lccdf(OpenCL)";
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
  auto one_m_erf
      = select(y_val < -37.5, 2.0,
               select(y_val < -5.0, 2.0 - erfc(-scaled_y),
                      select(y_val > 8.25, 0.0, 1.0 - erf(scaled_y))));
  auto lccdf_expr = colwise_sum(log(one_m_erf));
  auto y_deriv = -select(
      y_val > 8.25, INFTY,
      SQRT_TWO_OVER_SQRT_PI * elt_divide(exp(-square(scaled_y)), one_m_erf));

  matrix_cl<double> lccdf_cl;
  matrix_cl<double> y_deriv_cl;

  results(check_y_not_nan, lccdf_cl, y_deriv_cl)
      = expressions(y_not_nan_expr, lccdf_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv));

  T_partials_return lccdf = from_matrix_cl(lccdf_cl).sum() + LOG_HALF * N;

  operands_and_partials<decltype(y_col)> ops_partials(y_col);

  if (!is_constant<T_y_cl>::value) {
    ops_partials.edge1_.partials_ = std::move(y_deriv_cl);
  }
  return ops_partials.build(lccdf);
}

}  // namespace math
}  // namespace stan
#endif
#endif
