#ifndef STAN_MATH_OPENCL_GP_EXP_QUAD_COV_HPP
#define STAN_MATH_OPENCL_GP_EXP_QUAD_COV_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernels/gp_exp_quad_cov.hpp>
#include <stan/math/opencl/err/check_matching_dims.hpp>
#include <CL/cl.hpp>

namespace stan {
namespace math {
/**
 * Squared exponential kernel on the GPU.
 *
 * @param x input vector or matrix
 * @param sigma standard deviation
 * @param length_scale length scale
 *
 * @return Squared distance between elements of x.
 */
inline matrix_cl gp_exp_quad_cov(const matrix_cl& x, const double sigma,
                                 const double length_scale) {
  matrix_cl res(x.cols(), x.cols());
  try {
    opencl_kernels::gp_exp_quad_cov(cl::NDRange(x.cols(), x.cols()), x, res,
                                    sigma * sigma, -0.5 / square(length_scale),
                                    x.cols(), x.rows());
  } catch (const cl::Error& e) {
    check_opencl_error("gp_exp_quad_cov", e);
  }
  return res;
}

/**
 * Squared exponential kernel on the GPU.
 *
 * This function is for the cross covariance
 * matrix needed to compute the posterior predictive density.
 *
 * @param x first input vector or matrix
 * @param y second input vector or matrix
 * @param sigma standard deviation
 * @param length_scale length scale
 *
 * @return Squared distance between elements of x and y.
 */
inline matrix_cl gp_exp_quad_cov(const matrix_cl& x, const matrix_cl& y,
                                 const double sigma,
                                 const double length_scale) {
  check_size_match("gp_exp_quad_cov_cross", "x", x.rows(), "y", y.rows());
  matrix_cl res(x.cols(), y.cols());
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
