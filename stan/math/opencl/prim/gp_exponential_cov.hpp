#ifndef STAN_MATH_OPENCL_PRIM_GP_EXPONENTIAL_COV_HPP
#define STAN_MATH_OPENCL_PRIM_GP_EXPONENTIAL_COV_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernels/gp_exponential_cov.hpp>
#include <stan/math/opencl/err.hpp>
#include <CL/opencl.hpp>

namespace stan {
namespace math {
/** \ingroup opencl
 * Matern exponential kernel on the GPU.
 *
 * @tparam Matcl1 Type of the matrix
 * @tparam Arith2 Type of sigma
 * @tparam Arith3 Type of length_scale
 * @param x input vector or matrix
 * @param sigma standard deviation
 * @param length_scale length scale
 *
 * @return dot product covariance matrix that is positive semi-definite
 */
template <typename Matcl1, typename Arith2, typename Arith3,
          require_all_kernel_expressions_and_none_scalar_t<Matcl1>* = nullptr,
          require_all_arithmetic_t<Arith2, Arith3>* = nullptr>
inline matrix_cl<return_type_t<Matcl1, Arith2, Arith3>> gp_exponential_cov(
    const Matcl1& x, const Arith2 sigma, const Arith3 length_scale) {
  const auto& x_eval = x.eval();
  matrix_cl<return_type_t<Matcl1, Arith2, Arith3>> res(x.cols(), x.cols());
  int block_size = 16;
  int n_blocks = (x.cols() + block_size - 1) / block_size;
  int blocked_size = block_size * n_blocks;
  try {
    opencl_kernels::gp_exponential_cov(cl::NDRange(blocked_size, blocked_size),
                                       cl::NDRange(block_size, block_size),
                                       x_eval, res, sigma * sigma,
                                       -1.0 / length_scale, x.cols(), x.rows());
  } catch (const cl::Error& e) {
    check_opencl_error("gp_exponential_cov", e);
  }
  return res;
}

/** \ingroup opencl
 * Matern exponential kernel on the GPU.
 *
 * This function is for the cross covariance
 * matrix needed to compute the posterior predictive density.
 *
 * @tparam Matcl1 Type of the first matrix
 * @tparam Matcl2 Type of the second matrix
 * @tparam Arith3 Type of sigma
 * @tparam Arith4 Type of length scale
 * @param x first input vector or matrix
 * @param y second input vector or matrix
 * @param sigma standard deviation
 * @param length_scale length scale
 *
 * @return dot product covariance matrix
 */
template <typename Matcl1, typename Matcl2, typename Arith3, typename Arith4,
          require_all_kernel_expressions_and_none_scalar_t<Matcl1, Matcl2>* = nullptr,
          require_all_arithmetic_t<Arith3, Arith4>* = nullptr>
inline matrix_cl<return_type_t<Matcl1, Matcl2, Arith3, Arith4>> gp_exponential_cov(
    const Matcl1& x, const Matcl2& y, const Arith3 sigma, const Arith4 length_scale) {
  check_size_match("gp_exponential_cov_cross", "x", x.rows(), "y", y.rows());
  matrix_cl<return_type_t<Matcl1, Matcl2, Arith3, Arith4>> res(x.cols(), y.cols());
  const auto& x_eval = x.eval();
  const auto& y_eval = y.eval();
  int block_size = 16;
  int x_blocks = (x.cols() + block_size - 1) / block_size;
  int x_blocked_size = block_size * x_blocks;
  int y_blocks = (y.cols() + block_size - 1) / block_size;
  int y_blocked_size = block_size * y_blocks;
  try {
    opencl_kernels::gp_exponential_cov_cross(
        cl::NDRange(x_blocked_size, y_blocked_size),
        cl::NDRange(block_size, block_size), x_eval, y_eval, res, sigma * sigma,
        -1.0 / length_scale, x.cols(), y.cols(), x.rows());
  } catch (const cl::Error& e) {
    check_opencl_error("gp_exponential_cov_cross", e);
  }
  return res;
}

/** \ingroup opencl
 * Squared exponential kernel on the GPU.
 *
 * @tparam Matcl1 Type of the matrix
 * @tparam Arith2 Type of sigma
 * @tparam Matcl3 Type of length_scale
 * @param x input vector or matrix
 * @param sigma standard deviation
 * @param length_scale length scale
 *
 * @return Squared distance between elements of x.
 */
template <typename Matcl1, typename Arith2, typename Matcl3,
          require_all_kernel_expressions_and_none_scalar_t<Matcl1, Matcl3>* = nullptr,
          require_all_arithmetic_t<Arith2>* = nullptr>
inline matrix_cl<return_type_t<Matcl1, Arith2, Matcl3>> gp_exponential_cov(
    const Matcl1& x, const Arith2 sigma, const Matcl3& length_scale) {
  const auto& x_eval = elt_divide(x, rowwise_broadcast(length_scale)).eval();
  matrix_cl<return_type_t<Matcl1, Arith2, Matcl3>> res(x.cols(), x.cols());
  int block_size = 16;
  int n_blocks = (x.cols() + block_size - 1) / block_size;
  int blocked_size = block_size * n_blocks;
  try {
    opencl_kernels::gp_exponential_cov(cl::NDRange(blocked_size, blocked_size),
                                       cl::NDRange(block_size, block_size),
                                       x_eval, res, sigma * sigma, -1.0,
                                       x.cols(), x.rows());
  } catch (const cl::Error& e) {
    check_opencl_error("gp_exponential_cov", e);
  }
  return res;
}

/** \ingroup opencl
 * Squared exponential kernel on the GPU.
 *
 * This function is for the cross covariance
 * matrix needed to compute the posterior predictive density.
 *
 * @tparam Matcl1 Type of the first matrix
 * @tparam Matcl2 Type of the second matrix
 * @tparam Arith3 Type of sigma
 * @tparam Matcl4 Type of length scale
 * @param x first input vector or matrix
 * @param y second input vector or matrix
 * @param sigma standard deviation
 * @param length_scale length scale
 *
 * @return Squared distance between elements of x and y.
 */
template <
    typename Matcl1, typename Matcl2, typename Arith3, typename Matcl4,
    require_all_kernel_expressions_and_none_scalar_t<Matcl1, Matcl2, Matcl4>* = nullptr,
    require_all_arithmetic_t<Arith3>* = nullptr>
inline matrix_cl<return_type_t<Matcl1, Matcl2, Matcl3, Matcl4>> gp_exponential_cov(
    const Matcl1& x, const Matcl2& y, const Arith3 sigma, const Matcl4& length_scale) {
  check_size_match("gp_exponential_cov_cross", "x", x.rows(), "y", y.rows());
  matrix_cl<return_type_t<Matcl1, Matcl2, Arith3, Matcl4>> res(x.cols(), y.cols());
  const auto& x_eval = elt_divide(x, rowwise_broadcast(length_scale)).eval();
  const auto& y_eval = elt_divide(y, rowwise_broadcast(length_scale)).eval();
  int block_size = 16;
  int x_blocks = (x.cols() + block_size - 1) / block_size;
  int x_blocked_size = block_size * x_blocks;
  int y_blocks = (y.cols() + block_size - 1) / block_size;
  int y_blocked_size = block_size * y_blocks;
  try {
    opencl_kernels::gp_exponential_cov_cross(
        cl::NDRange(x_blocked_size, y_blocked_size),
        cl::NDRange(block_size, block_size), x_eval, y_eval, res, sigma * sigma,
        -1.0, x.cols(), y.cols(), x.rows());
  } catch (const cl::Error& e) {
    check_opencl_error("gp_exponential_cov_cross", e);
  }
  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
