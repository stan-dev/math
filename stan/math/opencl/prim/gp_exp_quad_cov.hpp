#ifndef STAN_MATH_OPENCL_PRIM_GP_EXP_QUAD_COV_HPP
#define STAN_MATH_OPENCL_PRIM_GP_EXP_QUAD_COV_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernels/gp_exp_quad_cov.hpp>
#include <stan/math/opencl/err.hpp>
#include <CL/opencl.hpp>

namespace stan {
namespace math {
/** \ingroup opencl
 * Squared exponential kernel on the GPU.
 *
 * @tparam T1 Type of the matrix
 * @tparam T2 Type of sigma
 * @tparam T3 Type of length_scale
 * @param x input vector or matrix
 * @param sigma standard deviation
 * @param length_scale length scale
 *
 * @return Squared distance between elements of x.
 */
template <typename T1, typename T2, typename T3,
          typename = require_all_arithmetic_t<T1, T2, T3>>
inline matrix_cl<return_type_t<T1, T2, T3>> gp_exp_quad_cov(
    const matrix_cl<T1>& x, const T2 sigma, const T3 length_scale) {
  matrix_cl<return_type_t<T1, T2, T3>> res(x.cols(), x.cols());
  try {
    opencl_kernels::gp_exp_quad_cov(cl::NDRange(x.cols(), x.cols()), x, res,
                                    sigma * sigma, -0.5 / square(length_scale),
                                    x.cols(), x.rows());
  } catch (const cl::Error& e) {
    check_opencl_error("gp_exp_quad_cov", e);
  }
  return res;
}

/** \ingroup opencl
 * Squared exponential kernel on the GPU.
 *
 * This function is for the cross covariance
 * matrix needed to compute the posterior predictive density.
 *
 * @tparam T1 Type of the first matrix
 * @tparam T2 Type of the second matrix
 * @tparam T3 Type of sigma
 * @tparam T4 Type of length scale
 * @param x first input vector or matrix
 * @param y second input vector or matrix
 * @param sigma standard deviation
 * @param length_scale length scale
 *
 * @return Squared distance between elements of x and y.
 */
template <typename T1, typename T2, typename T3, typename T4,
          typename = require_all_arithmetic_t<T1, T2, T3, T4>>
inline matrix_cl<return_type_t<T1, T2, T3, T4>> gp_exp_quad_cov(
    const matrix_cl<T1>& x, const matrix_cl<T2>& y, const T3 sigma,
    const T4 length_scale) {
  check_size_match("gp_exp_quad_cov_cross", "x", x.rows(), "y", y.rows());
  matrix_cl<return_type_t<T1, T2, T3, T4>> res(x.cols(), y.cols());
  try {
    opencl_kernels::gp_exp_quad_cov_cross(
        cl::NDRange(x.cols(), y.cols()), x, y, res, sigma * sigma,
        -0.5 / square(length_scale), x.cols(), y.cols(), x.rows());
  } catch (const cl::Error& e) {
    check_opencl_error("gp_exp_quad_cov_cross", e);
  }
  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
