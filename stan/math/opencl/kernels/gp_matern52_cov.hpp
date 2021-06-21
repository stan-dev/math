#ifndef STAN_MATH_OPENCL_KERNELS_gp_MATERN52_COV_HPP
#define STAN_MATH_OPENCL_KERNELS_gp_MATERN52_COV_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static const std::string gp_matern52_cov_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     * GPU part of calculation of Matern 5/2 kernel.
     *
     * @param[in] x input vector or matrix
     * @param[out] res squared distances between elements of x
     * @param sigma_sq squared standard deviation
     * @param root_5_inv_l sqrt(5.0) / length_scale
     * @param inv_l_sq_5_3 5.0 / 3.0 / square(length_scale)
     * @param size number of elements in x
     * @param element_size the number of doubles that make one element of x
     */
    __kernel void gp_matern52_cov(
        const __global double* x, __global double* res, const double sigma_sq,
        const double root_5_inv_l, const double inv_l_sq_5_3, const int size,
        const int element_size) {
      const int i = get_global_id(0);
      const int j = get_global_id(1);
      if (i < size && j < size) {
        if (i > j) {
          double sum = 0;
          for (int k = 0; k < element_size; k++) {
            double d = x[i * element_size + k] - x[j * element_size + k];
            sum += d * d;
          }
          double dist = sqrt(sum);
          double a = sigma_sq * (1.0 + root_5_inv_l * dist + inv_l_sq_5_3 * sum)
                     * exp(-root_5_inv_l * dist);
          res[j * size + i] = a;
          res[i * size + j] = a;
        } else if (i == j) {
          res[j * size + i] = sigma_sq;
        }
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/gp_matern52_cov.hpp gp_matern52_cov() \endlink
 */
const kernel_cl<in_buffer, out_buffer, double, double, double, int, int>
    gp_matern52_cov("gp_matern52_cov", {gp_matern52_cov_kernel_code});

// \cond
static const std::string gp_matern52_cov_cross_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     * GPU part of calculation of Matern 5/2 kernel.
     *
     * This function is for the cross covariance
     * matrix needed to compute the posterior predictive density.
     *
     * @param[in] x1 first input vector or matrix
     * @param[in] x2 second input vector or matrix
     * @param[out] res squared distances between elements of x
     * @param sigma_sq squared standard deviation
     * @param root_5_inv_l sqrt(5.0) / length_scale
     * @param inv_l_sq_5_3 5.0 / 3.0 / square(length_scale)
     * @param size1 number of elements in x1
     * @param size2 number of elements in x2
     * @param element_size the number of doubles that make one element of x and
     * y
     */
    __kernel void gp_matern52_cov_cross(
        const __global double* x1, const __global double* x2,
        __global double* res, const double sigma_sq, const double root_5_inv_l,
        const double inv_l_sq_5_3, const int size1, const int size2,
        const int element_size) {
      const int i = get_global_id(0);
      const int j = get_global_id(1);
      if (i < size1 && j < size2) {
        double sum = 0;
        for (int k = 0; k < element_size; k++) {
          double d = x1[i * element_size + k] - x2[j * element_size + k];
          sum += d * d;
        }
        double dist = sqrt(sum);
        res[j * size1 + i] = sigma_sq
                             * (1.0 + root_5_inv_l * dist + inv_l_sq_5_3 * sum)
                             * exp(-root_5_inv_l * dist);
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/gp_matern52_cov.hpp gp_matern52_cov_cross()
 * \endlink
 */
const kernel_cl<in_buffer, in_buffer, out_buffer, double, double, double, int,
                int, int>
    gp_matern52_cov_cross("gp_matern52_cov_cross",
                          {gp_matern52_cov_cross_kernel_code});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
