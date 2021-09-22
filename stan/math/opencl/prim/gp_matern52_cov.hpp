#ifndef STAN_MATH_OPENCL_PRIM_GP_MATERN52_COV_HPP
#define STAN_MATH_OPENCL_PRIM_GP_MATERN52_COV_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernels/gp_matern52_cov.hpp>
#include <stan/math/opencl/err.hpp>
#include <CL/opencl.hpp>

namespace stan {
namespace math {
/** \ingroup opencl
 * Matern 5/2 kernel on the GPU.
 *
 * @tparam T1 Type of the matrix
 * @tparam T2 Type of sigma
 * @tparam T3 Type of length_scale
 * @param x input vector or matrix
 * @param sigma standard deviation
 * @param length_scale length scale
 *
 * @return Matern 5/2 covariance matrix
 */
template <typename T1, typename T2, typename T3,
          require_all_kernel_expressions_and_none_scalar_t<T1>* = nullptr,
          require_all_arithmetic_t<T2, T3>* = nullptr>
inline matrix_cl<return_type_t<T1, T2, T3>> gp_matern52_cov(
    const T1& x, const T2 sigma, const T3 length_scale) {
  const auto& x_eval = x.eval();
  matrix_cl<return_type_t<T1, T2, T3>> res(x.cols(), x.cols());
  int block_size = 16;
  int n_blocks = (x.cols() + block_size - 1) / block_size;
  int blocked_size = block_size * n_blocks;
  try {
    opencl_kernels::gp_matern52_cov(
        cl::NDRange(blocked_size, blocked_size),
        cl::NDRange(block_size, block_size), x_eval, res, sigma * sigma,
        std::sqrt(5.0) / length_scale, 5.0 / (3.0 * square(length_scale)),
        x.cols(), x.rows());
  } catch (const cl::Error& e) {
    check_opencl_error("gp_matern52_cov", e);
  }
  return res;
}

/** \ingroup opencl
 * Matern 5/2 kernel on the GPU.
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
 * @return Matern 5/2 covariance matrix
 */
template <typename T1, typename T2, typename T3, typename T4,
          require_all_kernel_expressions_and_none_scalar_t<T1, T2>* = nullptr,
          require_all_arithmetic_t<T3, T4>* = nullptr>
inline matrix_cl<return_type_t<T1, T2, T3, T4>> gp_matern52_cov(
    const T1& x, const T2& y, const T3 sigma, const T4 length_scale) {
  check_size_match("gp_matern52_cov_cross", "x", x.rows(), "y", y.rows());
  matrix_cl<return_type_t<T1, T2, T3, T4>> res(x.cols(), y.cols());
  const auto& x_eval = x.eval();
  const auto& y_eval = y.eval();
  int block_size = 16;
  int x_blocks = (x.cols() + block_size - 1) / block_size;
  int x_blocked_size = block_size * x_blocks;
  int y_blocks = (y.cols() + block_size - 1) / block_size;
  int y_blocked_size = block_size * y_blocks;
  try {
    opencl_kernels::gp_matern52_cov_cross(
        cl::NDRange(x_blocked_size, y_blocked_size),
        cl::NDRange(block_size, block_size), x_eval, y_eval, res, sigma * sigma,
        std::sqrt(5.0) / length_scale, 5.0 / (3.0 * square(length_scale)),
        x.cols(), y.cols(), x.rows());
  } catch (const cl::Error& e) {
    check_opencl_error("gp_matern52_cov_cross", e);
  }
  return res;
}

/** \ingroup opencl
 * Matern 5/2 kernel on the GPU.
 *
 * @tparam T1 Type of the matrix
 * @tparam T2 Type of sigma
 * @tparam T3 Type of length_scale
 * @param x input vector or matrix
 * @param sigma standard deviation
 * @param length_scale length scale
 *
 * @return Matern 5/2 covariance matrix
 */
template <typename T1, typename T2, typename T3,
          require_all_kernel_expressions_and_none_scalar_t<T1, T3>* = nullptr,
          require_all_arithmetic_t<T2>* = nullptr>
inline matrix_cl<return_type_t<T1, T2, T3>> gp_matern52_cov(
    const T1& x, const T2 sigma, const T3 length_scale) {
  const auto& x_eval = elt_divide(x, rowwise_broadcast(length_scale)).eval();
  matrix_cl<return_type_t<T1, T2, T3>> res(x.cols(), x.cols());
  int block_size = 16;
  int n_blocks = (x.cols() + block_size - 1) / block_size;
  int blocked_size = block_size * n_blocks;
  try {
    opencl_kernels::gp_matern52_cov(cl::NDRange(blocked_size, blocked_size),
                                    cl::NDRange(block_size, block_size), x_eval,
                                    res, sigma * sigma, std::sqrt(5.0),
                                    5.0 / 3.0, x.cols(), x.rows());
  } catch (const cl::Error& e) {
    check_opencl_error("gp_matern52_cov", e);
  }
  return res;
}

/** \ingroup opencl
 * Matern 5/2 kernel on the GPU.
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
 * @return Matern 5/2 covariance matrix
 */
template <
    typename T1, typename T2, typename T3, typename T4,
    require_all_kernel_expressions_and_none_scalar_t<T1, T2, T4>* = nullptr,
    require_all_arithmetic_t<T3>* = nullptr>
inline matrix_cl<return_type_t<T1, T2, T3, T4>> gp_matern52_cov(
    const T1& x, const T2& y, const T3 sigma, const T4 length_scale) {
  check_size_match("gp_matern52_cov_cross", "x", x.rows(), "y", y.rows());
  matrix_cl<return_type_t<T1, T2, T3, T4>> res(x.cols(), y.cols());
  const auto& x_eval = elt_divide(x, rowwise_broadcast(length_scale)).eval();
  const auto& y_eval = elt_divide(y, rowwise_broadcast(length_scale)).eval();
  int block_size = 16;
  int x_blocks = (x.cols() + block_size - 1) / block_size;
  int x_blocked_size = block_size * x_blocks;
  int y_blocks = (y.cols() + block_size - 1) / block_size;
  int y_blocked_size = block_size * y_blocks;
  try {
    opencl_kernels::gp_matern52_cov_cross(
        cl::NDRange(x_blocked_size, y_blocked_size),
        cl::NDRange(block_size, block_size), x_eval, y_eval, res, sigma * sigma,
        std::sqrt(5.0), 5.0 / 3.0, x.cols(), y.cols(), x.rows());
  } catch (const cl::Error& e) {
    check_opencl_error("gp_matern52_cov_cross", e);
  }
  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
