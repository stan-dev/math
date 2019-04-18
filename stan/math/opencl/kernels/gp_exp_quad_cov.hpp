#ifndef STAN_MATH_OPENCL_KERNELS_GP_EXP_QUAD_COV_HPP
#define STAN_MATH_OPENCL_KERNELS_GP_EXP_QUAD_COV_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
const char* gp_exp_quad_cov_kernel_code = STRINGIFY(
// \endcond
/**
 *
 */
        __kernel void gp_exp_quad_cov(__global double* x, __global double* res, double sigma_sq, double neg_half_inv_l_sq, int size, int element_size) {
          const int i = get_global_id(0);
          const int j = get_global_id(1);
          if (i < size && j < (size - 1) && i > j) {
            double sum=0;
            for(int k=0;k<element_size;k++){
              double d = x[i * element_size + k] - x[j * element_size + k];
              sum += d*d;
            }
            double a = sigma_sq * exp(neg_half_inv_l_sq * sum);
            res[j * size + i] = a;
            res[i * size + j] = a;
          }
          else if (i == j) {
              res[j * size + i] = sigma_sq;
          }
        }
// \cond
);
// \endcond

// \cond
const char* gp_exp_quad_cov_2_kernel_code = STRINGIFY(
// \endcond
/**
 *
 */
        __kernel void gp_exp_quad_cov_2(__global double* x1, __global double* x2, __global double* res,
                                        double sigma_sq, double neg_half_inv_l_sq, int size1, int size2, int element_size) {
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

/**
 * See the docs for \link kernels/gp_exp_quad_cov.hpp add() \endlink
 */
const global_range_kernel<cl::Buffer, cl::Buffer, double, double, int, int>
        gp_exp_quad_cov("gp_exp_quad_cov", gp_exp_quad_cov_kernel_code,
                        {});


/**
 * See the docs for \link kernels/gp_exp_quad_cov.hpp add() \endlink
 */
const global_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, double, double, int, int, int>
        gp_exp_quad_cov_2("gp_exp_quad_cov_2", gp_exp_quad_cov_2_kernel_code,
                         {});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
