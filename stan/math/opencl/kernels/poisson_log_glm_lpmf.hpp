#ifndef STAN_MATH_OPENCL_KERNELS_POISSON_LOG_GLM_LPMF_HPP
#define STAN_MATH_OPENCL_KERNELS_POISSON_LOG_GLM_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static const char* poisson_log_glm_kernel_code = STRINGIFY(
    // \endcond
    /**
     * GPU implementation of Generalized Linear Model (GLM)
     * with Poisson distribution and log link function.
     *
     * Must be run with at least N threads and local size equal to LOCAL_SIZE_.
     * @param[out] theta_derivative_global derivative with respect to x * beta +
     * alpha
     * @param[out] theta_derivative_sum partially summed theta_derivative_global
     * (1 value per work group)
     * @param[out] logp_global partially summed part of log probability (1 value
     * per work group)
     * @param[in] y_global positive integer vector parameter
     * @param[in] x design matrix
     * @param[in] alpha intercept (in log odds)
     * @param[in] beta weight vector
     * @param N number of cases
     * @param M number of attributes
     * @param is_alpha_vector 0 or 1 - whether alpha is a vector (alternatively
     * it is a scalar)
     * @param need_logp1 interpreted as boolean - whether first part of
     * logp_global needs to be computed
     * @param need_logp2 interpreted as boolean - whether second part of
     * logp_global needs to be computed
     */
    __kernel void poisson_log_glm(
        __global double* theta_derivative_global,
        __global double* theta_derivative_sum, __global double* logp_global,
        const __global int* y_global, const __global double* x,
        const __global double* alpha, const __global double* beta, const int N,
        const int M, const int is_alpha_vector, const int need_logp1,
        const int need_logp2) {
      const int gid = get_global_id(0);
      const int lid = get_local_id(0);
      const int lsize = get_local_size(0);
      const int wg_id = get_group_id(0);

      __local double local_storage[LOCAL_SIZE_];
      double theta = 0;
      double theta_derivative = 0;
      double logp = 0;
      // Most calculations only happen for relevant data within next if.
      // Exceptions are reductions between threads that need barriers.
      if (gid < N) {
        for (int i = 0, j = 0; i < M; i++, j += N) {
          theta += x[j + gid] * beta[i];
        }

        theta += alpha[gid * is_alpha_vector];
        const double y = y_global[gid];
        const double exp_theta = exp(theta);
        theta_derivative = y - exp_theta;
        if (y < 0 || !isfinite(theta)) {
          // this signals that an exception must be raised
          theta_derivative = NAN;
        }
        if (need_logp1) {
          logp = -lgamma(y + 1);
        }
        if (need_logp2) {
          logp += y * theta - exp_theta;
        }
        theta_derivative_global[gid] = theta_derivative;
      }
      // Sum theta_derivative, calculated by different threads.
      // Since we can't sum between different work groups, we emit one number
      // per work group. These must be summed on CPU for final result.
      local_storage[lid] = theta_derivative;
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
        theta_derivative_sum[wg_id] = local_storage[0];
      }

      if (need_logp1 || need_logp2) {
        // Sum logp, calculated by different threads.
        barrier(CLK_LOCAL_MEM_FENCE);
        local_storage[lid] = logp;
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
          logp_global[wg_id] = local_storage[0];
        }
      }
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/subtract.hpp subtract() \endlink
 */
const kernel_cl<out_buffer, out_buffer, out_buffer, in_buffer, in_buffer,
                in_buffer, in_buffer, int, int, int, int, int>
    poisson_log_glm("poisson_log_glm", {poisson_log_glm_kernel_code},
                    {{"REDUCTION_STEP_SIZE", 4}, {"LOCAL_SIZE_", 64}});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
