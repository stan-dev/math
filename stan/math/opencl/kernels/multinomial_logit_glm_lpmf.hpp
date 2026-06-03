#ifndef STAN_MATH_OPENCL_KERNELS_MULTINOMIAL_LOGIT_GLM_LPMF_HPP
#define STAN_MATH_OPENCL_KERNELS_MULTINOMIAL_LOGIT_GLM_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static constexpr const char* multinomial_logit_glm_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     * GPU kernel for the Generalized Linear Model (GLM) with multinomial
     * distribution and softmax (logit) link function.
     *
     * Must be run with at least N_instances threads and local size LOCAL_SIZE_.
     * Each thread handles one instance n = gid.
     *
     * The kernel performs two passes over the K classes for instance n:
     *   1. find max(eta[n,:]) for numerical stability,
     *   2. accumulate sum_exp, S_n, and logp using shifted eta
     *      (eta[n,k] - max) to avoid catastrophic cancellation; skips
     *      skips y_nk=0 terms to implement the 0*log(0)=0 convention;
     *      if need_delta, stash exp(eta[n,k] - max) into delta_global.
     * A final loop normalizes delta (if need_delta) and subtracts
     * lgamma(y_nk+1) terms (if need_logp_gamma), reading only y_global
     * and delta_global, without re-reading x_beta_global or alpha_global.
     *
     * @param[out] logp_global partial logp sums, one per work group
     * @param[out] delta_global residual matrix N_instances x N_classes
     * (col-major)
     * @param[in] y_global outcome counts, N_instances x N_classes (col-major)
     * @param[in] x_beta_global product x*beta, N_instances x N_classes
     * (col-major)
     * @param[in] alpha_global intercepts: K values if is_alpha_vector, else
     *   N_instances x N_classes (col-major)
     * @param N_instances number of instances
     * @param N_classes number of outcome classes
     * @param is_alpha_vector 1 if alpha is shared 1xK row, 0 if NxK
     * @param need_delta 1 if delta_global should be computed and written
     * @param need_logp_gamma 1 if lgamma terms should be included in logp
     */
    __kernel void multinomial_logit_glm(
        __global double* logp_global, __global double* delta_global,
        const __global int* y_global, const __global double* x_beta_global,
        const __global double* alpha_global, const int N_instances,
        const int N_classes, const int is_alpha_vector, const int need_delta,
        const int need_logp_gamma) {
      const int gid = get_global_id(0);
      const int lid = get_local_id(0);
      const int lsize = get_local_size(0);
      const int wg_id = get_group_id(0);

      __local double local_storage[LOCAL_SIZE_];

      double logp = 0;
      if (gid < N_instances) {
        // Pass 1: row-wise max of eta for numerical stability.
        double eta_max = -DBL_MAX;
        for (int k = 0; k < N_classes; k++) {
          int nk = k * N_instances + gid;
          int alpha_idx = is_alpha_vector * k + !is_alpha_vector * nk;
          double eta_k = x_beta_global[nk] + alpha_global[alpha_idx];
          eta_max = fmax(eta_k, eta_max);
        }

        // Pass 2: sum_exp, S_n, logp; if need_delta stash exp_k in
        // delta_global.
        double sum_exp = 0;
        int S_n = 0;
        for (int k = 0; k < N_classes; k++) {
          int nk = k * N_instances + gid;
          int alpha_idx = is_alpha_vector * k + !is_alpha_vector * nk;
          double shifted_eta_k
              = x_beta_global[nk] + alpha_global[alpha_idx] - eta_max;
          double exp_k = exp(shifted_eta_k);
          sum_exp += exp_k;
          int y_nk = y_global[nk];
          S_n += y_nk;
          logp += (y_nk != 0) ? y_nk * shifted_eta_k : 0;
          if (need_delta) {
            delta_global[nk] = exp_k;
          }
        }
        logp -= S_n * log(sum_exp);

        if (need_logp_gamma) {
          logp += lgamma(S_n + 1.0);
        }

        // Normalize delta and/or subtract lgamma(y_nk+1) in one pass.
        if (need_delta || need_logp_gamma) {
          for (int k = 0; k < N_classes; k++) {
            int nk = k * N_instances + gid;
            int y_nk = y_global[nk];
            logp -= need_logp_gamma * lgamma(y_nk + 1.0);
            if (need_delta) {
              delta_global[nk] = y_nk - S_n * delta_global[nk] / sum_exp;
            }
          }
        }
      }

      // Work-group reduction of logp.
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
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/multinomial_logit_glm_lpmf.hpp
 * multinomial_logit_glm() \endlink
 */
const kernel_cl<out_buffer, out_buffer, in_buffer, in_buffer, in_buffer, int,
                int, int, int, int>
    multinomial_logit_glm("multinomial_logit_glm",
                          {multinomial_logit_glm_kernel_code},
                          {{"REDUCTION_STEP_SIZE", 4}, {"LOCAL_SIZE_", 64}});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
