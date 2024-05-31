#ifndef STAN_MATH_OPENCL_KERNELS_CATEGORICAL_LOGIT_GLM_LPMF_HPP
#define STAN_MATH_OPENCL_KERNELS_CATEGORICAL_LOGIT_GLM_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static constexpr const char* categorical_logit_glm_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     * GPU implementation of Generalized Linear Model (GLM)
     * with categorical distribution and logit (softmax) link function.
     *
     * Must be run with at least N threads and local size equal to LOCAL_SIZE_.
     * @param[out] logp_global partially summed log probability (1 value per
     * work group)
     * @param[out] exp_lin_global exponentiation of sum of alpha and matrix
     * product of x and beta
     * @param[out] inv_sum_exp_lin_global inverse of rowwise sum of \c
     * exp_lin_global
     * @param[out] neg_softmax_lin_global negated softmax of sum of alpha and
     * matrix product of x and beta
     * @param[out] alpha_derivative derivative wrt alpha
     * @param[in] y_global a scalar or vector of classes.
     * @param[in] x_beta_global product of design matrix and weight matrix
     * @param[in] alpha_global intercept (in log odds)
     * @param N_instances number of instances
     * @param N_attributes number of attributes
     * @param N_classes number of classes
     * @param is_y_vector 0 or 1 - whether y is a vector (alternatively it is a
     * scalar)
     * @param need_alpha_derivative interpreted as boolean - whether
     * alpha_derivative needs to be computed
     * @param need_neg_softmax_lin_global interpreted as boolean - whether
     * neg_softmax_lin_global needs to be computed
     */
    __kernel void categorical_logit_glm(
        __global double* logp_global, __global double* exp_lin_global,
        __global double* inv_sum_exp_lin_global,
        __global double* neg_softmax_lin_global,
        __global double* alpha_derivative, const __global int* y_global,
        const __global double* x_beta_global,
        const __global double* alpha_global, const int N_instances,
        const int N_attributes, const int N_classes, const int is_y_vector,
        const int need_alpha_derivative,
        const int need_neg_softmax_lin_global) {
      const int gid = get_global_id(0);
      const int lid = get_local_id(0);
      const int lsize = get_local_size(0);
      const int wg_id = get_group_id(0);
      const int ngroups = get_num_groups(0);

      __local double local_storage[LOCAL_SIZE_];

      double logp = 0;
      double inv_sum_exp_lin;
      int class_idx = -1;
      // Most calculations only happen for relevant data within next if.
      // Exceptions are reductions between threads that need barriers.
      if (gid < N_instances) {
        double lin_max = -INFINITY;
        for (int i = 0; i < N_classes; i++) {
          double lin = x_beta_global[i * N_instances + gid] + alpha_global[i];
          if (lin > lin_max) {
            lin_max = lin;
          }
        }
        double alpha = alpha_global[gid];
        double sum_exp_lin = 0;
        for (int i = 0; i < N_classes; i++) {
          double lin = x_beta_global[i * N_instances + gid] + alpha_global[i];
          double exp_lin = exp(lin - lin_max);
          sum_exp_lin += exp_lin;
          exp_lin_global[i * N_instances + gid] = exp_lin;
        }
        inv_sum_exp_lin = 1 / sum_exp_lin;
        inv_sum_exp_lin_global[gid] = inv_sum_exp_lin;

        class_idx = y_global[gid * is_y_vector] - 1;
        if (class_idx < 0 || class_idx > N_classes) {
          logp = NAN;
        } else {
          logp = log(inv_sum_exp_lin) - lin_max
                 + x_beta_global[class_idx * N_instances + gid]
                 + alpha_global[class_idx];
        }
      }
      barrier(CLK_GLOBAL_MEM_FENCE);
      double neg_softmax_lin_sum = 0;
      if (need_alpha_derivative || need_neg_softmax_lin_global) {
        for (int i = 0; i < N_classes; i++) {
          double neg_softmax_lin = 0;
          if (gid < N_instances) {
            int idx = i * N_instances + gid;
            neg_softmax_lin = -exp_lin_global[idx] * inv_sum_exp_lin;
            if (need_neg_softmax_lin_global) {
              neg_softmax_lin_global[idx] = neg_softmax_lin;
            }
          }
          if (need_alpha_derivative) {
            local_storage[lid] = neg_softmax_lin + (class_idx == i);
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
              alpha_derivative[i + wg_id * N_classes] = local_storage[0];
            }
            barrier(CLK_LOCAL_MEM_FENCE);
          }
        }
      }
      // Sum logp, calculated by different threads.
      // Since we can't sum between different work groups, we emit one number
      // per work group. These must be summed on CPU for final result.
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
 * See the docs for \link kernels/categorical_logit_glm_lpmf.hpp
 * categorical_logit_glm() \endlink
 */
const kernel_cl<out_buffer, out_buffer, out_buffer, out_buffer, out_buffer,
                in_buffer, in_buffer, in_buffer, int, int, int, int, int, int>
    categorical_logit_glm("categorical_logit_glm",
                          {categorical_logit_glm_kernel_code},
                          {{"REDUCTION_STEP_SIZE", 4}, {"LOCAL_SIZE_", 64}});

// \cond
static const std::string categorical_logit_glm_beta_derivative_kernel_code
    = STRINGIFY(
        // \endcond
        /** \ingroup opencl_kernels
         * Calculates derivative wrt beta.
         *
         * Must be run with global size of local_size*N_attributes.
         * @param[in,out] beta_derivative derivative wrt beta
         * @param temp temporary workspace of size global_size*N_classes
         * @param[in] y a scalar or vector of classes
         * @param[in] x design matrix
         * @param N_instances number of instances
         * @param N_attributes number of attributes
         * @param N_classes number of classes
         * @param is_y_vector 0 or 1 - whether y is a vector (alternatively it
         * is a scalar)
         */
        __kernel void categorical_logit_glm_beta_derivative(
            __global double* beta_derivative, __global double* temp,
            const __global int* y, const __global double* x,
            const int N_instances, const int N_attributes, const int N_classes,
            const int is_y_vector) {
          const int gid = get_global_id(0);
          const int lid = get_local_id(0);
          const int lsize = get_local_size(0);
          const int wg_id = get_group_id(0);

          for (int i = 0; i < N_classes; i++) {
            temp[gid * N_classes + i] = 0;
          }
          for (int i = lid; i < N_instances; i += lsize) {
            int pos = y[i * is_y_vector] - 1;
            temp[gid * N_classes + pos] += x[wg_id * N_instances + i];
          }
          barrier(CLK_GLOBAL_MEM_FENCE);
          for (int i = lid; i < N_classes; i += lsize) {
            double res = 0;
            for (int j = 0; j < lsize; j++) {
              res += temp[(wg_id * lsize + j) * N_classes + i];
            }
            beta_derivative[i * N_attributes + wg_id] += res;
          }
        }
        // \cond
    );  // NOLINT
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/categorical_logit_glm_lpmf.hpp
 * categorical_logit_glm_beta_derivative() \endlink
 */
const kernel_cl<in_out_buffer, in_out_buffer, in_buffer, in_buffer, int, int,
                int, int>
    categorical_logit_glm_beta_derivative(
        "categorical_logit_glm_beta_derivative",
        {categorical_logit_glm_beta_derivative_kernel_code});

}  // namespace opencl_kernels

}  // namespace math
}  // namespace stan
#endif
#endif
