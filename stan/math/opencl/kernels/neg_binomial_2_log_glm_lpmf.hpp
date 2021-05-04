#ifndef STAN_MATH_OPENCL_KERNELS_NEG_BINOMIAL_2_LOG_GLM_LPMF_HPP
#define STAN_MATH_OPENCL_KERNELS_NEG_BINOMIAL_2_LOG_GLM_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/kernels/device_functions/digamma.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static const char* neg_binomial_2_log_glm_kernel_code = STRINGIFY(
    // \endcond

    /** \ingroup opencl_kernels
     * GPU implementation of Generalized Linear Model (GLM)
     * with Negative-Binomial-2 distribution and log link function.
     *
     * Must be run with at least N threads and local size equal to LOCAL_SIZE_.
     * @param[out] logp_global partially summed log probability (1 value per
     * work group)
     * @param[out] theta_derivative_global intermediate variable used in the
     * model
     * @param[out] theta_derivative_sum partially summed theta_derivative_global
     * (1 value per work group)
     * @param[out] phi_derivative_global derivative with respect to phi
     * @param[in] y_global failures count vector parameter
     * @param[in] x design matrix
     * @param[in] alpha intercept (in log odds)
     * @param[in] beta weight vector
     * @param[in] phi_global (vector of) precision parameter(s)
     * @param N number of cases
     * @param M number of attributes
     * @param is_y_vector 0 or 1 - whether y is a vector (alternatively
     * it is a scalar)
     * @param is_alpha_vector 0 or 1 - whether alpha is a vector (alternatively
     * it is a scalar)
     * @param is_phi_vector 0 or 1 - whether phi is a vector (alternatively it
     * is a scalar)
     * @param need_theta_derivative whether theta_derivative needs to be
     * computed
     * @param need_theta_derivative_sum whether theta_derivative_sum needs to be
     * computed
     * @param need_phi_derivative whether phi_derivative needs to be computed
     * @param need_phi_derivative_sum whether phi_derivative_sum needs to be
     * computed
     * @param need_logp1 interpreted as boolean - whether first part logp_global
     * needs to be computed
     * @param need_logp2 interpreted as boolean - whether second part
     * logp_global needs to be computed
     * @param need_logp3 interpreted as boolean - whether third part logp_global
     * needs to be computed
     * @param need_logp4 interpreted as boolean - whether fourth part
     * logp_global needs to be computed
     */
    __kernel void neg_binomial_2_log_glm(
        __global double* logp_global, __global double* theta_derivative_global,
        __global double* theta_derivative_sum,
        __global double* phi_derivative_global, const __global int* y_global,
        const __global double* x, const __global double* alpha,
        const __global double* beta, const __global double* phi_global,
        const int N, const int M, const int is_y_vector,
        const int is_alpha_vector, const int is_phi_vector,
        const int need_theta_derivative, const int need_theta_derivative_sum,
        const int need_phi_derivative, const int need_phi_derivative_sum,
        const int need_logp1, const int need_logp2, const int need_logp3,
        const int need_logp4) {
      const int gid = get_global_id(0);
      const int lid = get_local_id(0);
      const int lsize = get_local_size(0);
      const int wgid = get_group_id(0);

      __local double res_loc[LOCAL_SIZE_];
      double logp = 0;
      double phi_derivative = 0;
      double theta_derivative = 0;

      // Most calculations only happen for relevant data within next if.
      // Exceptions are reductions between threads that need barriers.
      if (gid < N) {
        double theta = 0;
        for (int i = 0, j = 0; i < M; i++, j += N) {
          theta += x[j + gid] * beta[i];
        }
        double phi = phi_global[gid * is_phi_vector];
        double y = y_global[gid * is_y_vector];
        if (!isfinite(theta) || y < 0 || !isfinite(phi)) {
          logp = NAN;
        }
        theta += alpha[gid * is_alpha_vector];
        double log_phi = log(phi);
        double logsumexp_theta_logphi;
        if (theta > log_phi) {
          logsumexp_theta_logphi = theta + log1p(exp(log_phi - theta));
        } else {
          logsumexp_theta_logphi = log_phi + log1p(exp(theta - log_phi));
        }
        double y_plus_phi = y + phi;
        if (need_logp1) {
          logp -= lgamma(y + 1);
        }
        if (need_logp2) {
          logp -= lgamma(phi);
          if (phi != 0) {
            logp += phi * log(phi);
          }
        }
        logp -= y_plus_phi * logsumexp_theta_logphi;
        if (need_logp3) {
          logp += y * theta;
        }
        if (need_logp4) {
          logp += lgamma(y_plus_phi);
        }
        double theta_exp = exp(theta);
        theta_derivative = y - theta_exp * y_plus_phi / (theta_exp + phi);
        if (need_theta_derivative) {
          theta_derivative_global[gid] = theta_derivative;
        }
        if (need_phi_derivative) {
          phi_derivative = 1 - y_plus_phi / (theta_exp + phi) + log_phi
                           - logsumexp_theta_logphi + digamma(y_plus_phi)
                           - digamma(phi);
          if (!need_phi_derivative_sum) {
            phi_derivative_global[gid] = phi_derivative;
          }
        }
      }

      // Sum logp, calculated by different threads.
      // Since we can't sum between different work groups, we emit one number
      // per work group. These must be summed on CPU for final result.
      res_loc[lid] = logp;
      barrier(CLK_LOCAL_MEM_FENCE);
      for (int step = lsize / REDUCTION_STEP_SIZE; step > 0;
           step /= REDUCTION_STEP_SIZE) {
        if (lid < step) {
          for (int i = 1; i < REDUCTION_STEP_SIZE; i++) {
            res_loc[lid] += res_loc[lid + step * i];
          }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      if (lid == 0) {
        logp_global[wgid] = res_loc[0];
      }

      if (need_theta_derivative_sum) {
        // Sum theta_derivative, calculated by different threads.
        barrier(CLK_LOCAL_MEM_FENCE);
        res_loc[lid] = theta_derivative;
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int step = lsize / REDUCTION_STEP_SIZE; step > 0;
             step /= REDUCTION_STEP_SIZE) {
          if (lid < step) {
            for (int i = 1; i < REDUCTION_STEP_SIZE; i++) {
              res_loc[lid] += res_loc[lid + step * i];
            }
          }
          barrier(CLK_LOCAL_MEM_FENCE);
        }
        if (lid == 0) {
          theta_derivative_sum[wgid] = res_loc[0];
        }
      }

      if (need_phi_derivative_sum) {
        // Sum phi_derivative, calculated by different threads.
        barrier(CLK_LOCAL_MEM_FENCE);
        res_loc[lid] = phi_derivative;
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int step = lsize / REDUCTION_STEP_SIZE; step > 0;
             step /= REDUCTION_STEP_SIZE) {
          if (lid < step) {
            for (int i = 1; i < REDUCTION_STEP_SIZE; i++) {
              res_loc[lid] += res_loc[lid + step * i];
            }
          }
          barrier(CLK_LOCAL_MEM_FENCE);
        }
        if (lid == 0) {
          phi_derivative_global[wgid] = res_loc[0];
        }
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/neg_binomial_2_log_glm_lpmf.hpp
 * neg_binomial_2_log_glm_lpmf() \endlink
 */
const kernel_cl<out_buffer, out_buffer, out_buffer, out_buffer, in_buffer,
                in_buffer, in_buffer, in_buffer, in_buffer, int, int, int, int,
                int, int, int, int, int, int, int, int, int>
    neg_binomial_2_log_glm("neg_binomial_2_log_glm",
                           {digamma_device_function,
                            neg_binomial_2_log_glm_kernel_code},
                           {{"REDUCTION_STEP_SIZE", 4}, {"LOCAL_SIZE_", 64}});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
