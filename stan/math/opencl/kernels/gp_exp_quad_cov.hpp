#ifndef STAN_MATH_OPENCL_KERNELS_GP_EXP_QUAD_COV_HPP
#define STAN_MATH_OPENCL_KERNELS_GP_EXP_QUAD_COV_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static const std::string gp_exp_quad_cov_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     * GPU part of calculation of squared exponential kernel.
     *
     * @param[in] x input vector or matrix
     * @param[out] res squared distances between elements of x
     * @param sigma_sq squared standard deviation
     * @param neg_half_inv_l_sq -0.5 / square(length scale)
     * @param size number of elements in x
     * @param element_size the number of doubles that make one element of x
     */
    __kernel void gp_exp_quad_cov(const __global double* x,
                                  __global double* res, const double sigma_sq,
                                  const double neg_half_inv_l_sq,
                                  const int size, const int element_size) {
      const int i = get_global_id(0);
      const int j = get_global_id(1);
      if (i < size && j < (size - 1) && i > j) {
        double sum = 0;
        for (int k = 0; k < element_size; k++) {
          double d = x[i * element_size + k] - x[j * element_size + k];
          sum += d * d;
        }
        double a = sigma_sq * exp(neg_half_inv_l_sq * sum);
        res[j * size + i] = a;
        res[i * size + j] = a;
      } else if (i == j) {
        res[j * size + i] = sigma_sq;
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/gp_exp_quad_cov.hpp gp_exp_quad_cov() \endlink
 */
const kernel_cl<in_buffer, out_buffer, double, double, int, int>
    gp_exp_quad_cov("gp_exp_quad_cov", {gp_exp_quad_cov_kernel_code});

// \cond
static const std::string gp_exp_quad_cov_cross_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     * GPU part of calculation of squared exponential kernel.
     *
     * This function is for the cross covariance
     * matrix needed to compute the posterior predictive density.
     *
     * @param[in] x1 first input vector or matrix
     * @param[in] x2 second input vector or matrix
     * @param[out] res squared distances between elements of x
     * @param sigma_sq squared standard deviation
     * @param neg_half_inv_l_sq -0.5 / square(length scale)
     * @param size1 number of elements in x1
     * @param size2 number of elements in x2
     * @param element_size the number of doubles that make one element of x and
     * y
     */
    __kernel void gp_exp_quad_cov_cross(
        const __global double* x1, const __global double* x2,
        __global double* res, const double sigma_sq,
        const double neg_half_inv_l_sq, const int size1, const int size2,
        const int element_size) {
      const int i = get_global_id(0);
      const int j = get_global_id(1);
      if (i < size1 && j < size2) {
        double sum = 0;
        for (int k = 0; k < element_size; k++) {
          double d = x1[i * element_size + k] - x2[j * element_size + k];
          sum += d * d;
        }
        res[j * size1 + i] = sigma_sq * exp(neg_half_inv_l_sq * sum);
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/gp_exp_quad_cov.hpp gp_exp_quad_cov_cross()
 * \endlink
 */
const kernel_cl<in_buffer, in_buffer, out_buffer, double, double, int, int, int>
    gp_exp_quad_cov_cross("gp_exp_quad_cov_cross",
                          {gp_exp_quad_cov_cross_kernel_code});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
