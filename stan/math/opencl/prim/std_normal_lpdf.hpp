#ifndef STAN_MATH_OPENCL_PRIM_STD_NORMAL_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_STD_NORMAL_LPDF_HPP
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
 * The log of the normal density for the specified scalar(s) given
 * a location of 0 and a scale of 1. y can be either
 * a scalar or a vector.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each observation.
 *
 * @tparam T_y_cl type of scalar
 * @param y Sequence of scalars.
 * @return The log of the product of the densities.
 * @throw std::domain_error if any scalar is nan.
 */
template <bool propto, typename T_y_cl,
          require_all_prim_or_rev_kernel_expression_t<T_y_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl>* = nullptr>
inline return_type_t<T_y_cl> std_normal_lpdf(const T_y_cl& y) {
  static constexpr const char* function = "std_normal_lpdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl>;
  using std::isfinite;
  using std::isnan;

  const size_t N = math::size(y);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl>::value) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& y_val = value_of(y_col);

  auto check_y_not_nan
      = check_cl(function, "Random variable", y_val, "not NaN");
  auto y_not_nan = !isnan(y_val);

  auto logp_expr = colwise_sum(elt_multiply(y_val, y_val));

  auto y_deriv = -y_val;

  matrix_cl<double> logp_cl;
  matrix_cl<double> y_deriv_cl;

  results(check_y_not_nan, logp_cl, y_deriv_cl) = expressions(
      y_not_nan, logp_expr, calc_if<!is_constant<T_y_cl>::value>(y_deriv));

  T_partials_return logp = sum(from_matrix_cl(logp_cl)) * -0.5;

  if (include_summand<propto>::value) {
    logp += NEG_LOG_SQRT_TWO_PI * N;
  }

  auto ops_partials = make_partials_propagator(y_col);

  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
