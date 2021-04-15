#ifndef STAN_MATH_OPENCL_PRIM_GP_DOT_PROD_COV_HPP
#define STAN_MATH_OPENCL_PRIM_GP_DOT_PROD_COV_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/value_of.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Dot product kernel on the GPU.
 *
 * @tparam T1 Type of the matrix
 * @tparam T2 Type of sigma
 * @param x input matrix
 * @param sigma standard deviation
 * @param length_scale length scale
 *
 * @return dot product covariance matrix that is positive semi-definite
 */
template <typename T_x, typename T_sigma,
          require_all_prim_or_rev_kernel_expression_t<T_x>* = nullptr,
          require_stan_scalar_t<T_sigma>* = nullptr>
inline auto gp_dot_prod_cov(const T_x& x, const T_sigma sigma) {
  const char* fun = "gp_dot_prod_cov(OpenCL)";
  check_nonnegative(fun, "sigma", sigma);
  check_finite(fun, "sigma", sigma);
  const auto& x_val = value_of(x);
  check_cl(fun, "x", x_val, "not NaN") = !isnan(x_val);
  return add(square(sigma), transpose(x) * x);
}

/** \ingroup opencl
 * Dot product kernel on the GPU.
 *
 * @tparam T1 Type of the matrix
 * @tparam T2 Type of sigma
 * @param x input matrix
 * @param sigma standard deviation
 * @param length_scale length scale
 *
 * @return dot product covariance matrix
 */
template <typename T_x, typename T_y, typename T_sigma,
          require_all_prim_or_rev_kernel_expression_t<T_x, T_y>* = nullptr,
          require_stan_scalar_t<T_sigma>* = nullptr>
inline auto gp_dot_prod_cov(const T_x& x, const T_y& y, const T_sigma sigma) {
  const char* fun = "gp_dot_prod_cov(OpenCL)";
  check_nonnegative(fun, "sigma", sigma);
  check_finite(fun, "sigma", sigma);
  const auto& x_val = value_of(x);
  const auto& y_val = value_of(y);
  check_cl(fun, "x", x_val, "not NaN") = !isnan(x_val);
  check_cl(fun, "y", y_val, "not NaN") = !isnan(y_val);
  return add(square(sigma), transpose(x) * y);
}
}  // namespace math
}  // namespace stan
#endif
#endif
