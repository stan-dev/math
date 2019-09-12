#ifndef STAN_MATH_OPENCL_KERNELS_NORMAL_ID_GLM_LPDF_HPP
#define STAN_MATH_OPENCL_KERNELS_NORMAL_ID_GLM_LPDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static const char* normal_id_glm_kernel_code = STRINGIFY(
    // \endcond
    /**
     * GPU implementation of Generalized Linear Model (GLM)
     * with Normal distribution and identity link function.
     *
     * Must be run with at least N threads and local size equal to LOCAL_SIZE_.
     * @param[in] y vector parameter
     * @param[in] x design matrix
     * @param[in] alpha intercept (in log odds)
     * @param[in] beta weight vector
     * @param[in] sigma_global (Sequence of) scale parameters for the normal
     * @param[out] mu_derivative_global intermediate variable used in the model
     * @param[out] mu_derivative_sum partially summed mu_derivative_global (1
     * value per work group)
     * @param[out] y_scaled_sq_sum y, scaled by expected mean and given
     * variance, squared and partially summed (1 value per work group).
     * @param[out] sigma_derivative derivative with respect to sigma
     * @param[out] log_sigma_sum partially summed logarithm of sigma (1 value
     * per work group)
     * @param N number of cases
     * @param M number of attributes
     * @param is_alpha_vector 0 or 1 - whether alpha is a vector (alternatively
     * it is a scalar)
     * @param is_sigma_vector 0 or 1 - whether sigma is a vector (alternatively
     * it is a scalar)
     * @param need_mu_derivative interpreted as boolean - whether mu_derivative
     * needs to be computed
     * @param need_mu_derivative_sum interpreted as boolean - whether
     * mu_derivative_sum needs to be computed
     * @param need_sigma_derivative interpreted as boolean - whether
     * sigma_derivative needs to be computed
     * @param need_log_sigma_sum interpreted as boolean - whether log_sigma_sum
     * needs to be computed
     */
    __kernel void normal_id_glm(
        __global double* mu_derivative_global,
        __global double* mu_derivative_sum, __global double* y_scaled_sq_sum,
        __global double* sigma_derivative, __global double* log_sigma_sum,
        const __global double* y, const __global double* x,
        const __global double* alpha, const __global double* beta,
        const __global double* sigma_global, const int N, const int M,
        const int is_alpha_vector, const int is_sigma_vector,
        const int need_mu_derivative, const int need_mu_derivative_sum,
        const int need_sigma_derivative, const int need_log_sigma_sum) {
      const int gid = get_global_id(0);
      const int lid = get_local_id(0);
      const int lsize = get_local_size(0);
      const int wg_id = get_group_id(0);

      __local double local_storage[LOCAL_SIZE_];

      double y_scaled_sq = 0;
      double log_sigma = 0;
      double mu_derivative = 0;
      // Most calculations only happen for relevant data within next if.
      // Exceptions are reductions between threads that need barriers.
      if (gid < N) {
        double y_scaled = 0;
        for (int i = 0, j = 0; i < M; i++, j += N) {
          y_scaled += x[j + gid] * beta[i];
        }
        double sigma = sigma_global[gid * is_sigma_vector];
        double inv_sigma = 1 / sigma;
        y_scaled
            = (y[gid] - y_scaled - alpha[gid * is_alpha_vector]) * inv_sigma;
        mu_derivative = inv_sigma * y_scaled;
        if (need_mu_derivative) {
          mu_derivative_global[gid] = mu_derivative;
        }
        y_scaled_sq = y_scaled * y_scaled;
        if (need_sigma_derivative) {
          sigma_derivative[gid] = (y_scaled_sq - 1) * inv_sigma;
        }
        if (need_log_sigma_sum) {
          log_sigma = log(sigma);
        }
      }

      // Sum y_scaled_sq, calculated by different threads.
      // Since we can't sum between different work groups, we emit one number
      // per work group. These must be summed on CPU for final result.
      local_storage[lid] = y_scaled_sq;
      barrier(CLK_LOCAL_MEM_FENCE);
      for (int step = lsize / REDUCTION_STEP_SIZE; step > 0;
           step /= REDUCTION_STEP_SIZE) {
        if (lid < step) {
          for (int i = 1; i < REDUCTION_STEP_SIZE; i++) {
            local_storage[lid] += local_storage[lid + step * i];
          }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      if (lid == 0) {
        y_scaled_sq_sum[wg_id] = local_storage[0];
      }

      if (need_mu_derivative_sum) {
        // Sum mu_derivative, calculated by different threads.
        barrier(CLK_LOCAL_MEM_FENCE);
        local_storage[lid] = mu_derivative;
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int step = lsize / REDUCTION_STEP_SIZE; step > 0;
             step /= REDUCTION_STEP_SIZE) {
          if (lid < step) {
            for (int i = 1; i < REDUCTION_STEP_SIZE; i++) {
              local_storage[lid] += local_storage[lid + step * i];
            }
          }
          barrier(CLK_LOCAL_MEM_FENCE);
        }
        if (lid == 0) {
          mu_derivative_sum[wg_id] = local_storage[0];
        }
      }

      if (need_log_sigma_sum) {
        // Sum log_sigma, calculated by different threads.
        barrier(CLK_LOCAL_MEM_FENCE);
        local_storage[lid] = log_sigma;
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int step = lsize / REDUCTION_STEP_SIZE; step > 0;
             step /= REDUCTION_STEP_SIZE) {
          if (lid < step) {
            for (int i = 1; i < REDUCTION_STEP_SIZE; i++) {
              local_storage[lid] += local_storage[lid + step * i];
            }
          }
          barrier(CLK_LOCAL_MEM_FENCE);
        }
        if (lid == 0) {
          log_sigma_sum[wg_id] = local_storage[0];
        }
      }
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/normal_id_glm_lpdf.hpp
 * normal_id_glm() \endlink
 */
const kernel_cl<out_buffer, out_buffer, out_buffer, out_buffer, out_buffer,
                in_buffer, in_buffer, in_buffer, in_buffer, in_buffer, int, int,
                int, int, int, int, int, int>
    normal_id_glm("normal_id_glm", {normal_id_glm_kernel_code},
                  {{"REDUCTION_STEP_SIZE", 4}, {"LOCAL_SIZE_", 64}});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
