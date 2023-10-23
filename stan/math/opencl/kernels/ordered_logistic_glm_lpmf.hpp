#ifndef STAN_MATH_OPENCL_KERNELS_ORDERED_LOGISTIC_GLM_LPMF_HPP
#define STAN_MATH_OPENCL_KERNELS_ORDERED_LOGISTIC_GLM_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/kernels/device_functions/log1m_exp.hpp>
#include <stan/math/opencl/kernels/device_functions/log1p_exp.hpp>
#include <stan/math/opencl/kernels/device_functions/inv_logit.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static const char* ordered_logistic_glm_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     * GPU implementation of ordinal regression Generalized Linear Model (GLM).
     *
     * Must be run with at least N_instances threads and local size equal to
     * LOCAL_SIZE_.
     * @param[out] location_sum partially summed location (1 value per work
     * group)
     * @param[out] logp_global partially summed log probability (1 value per
     * work group)
     * @param[out] location_derivative derivative wrt location
     * @param[out] cuts_derivative partially summed derivative wrt cuts (1
     * column per work group)
     * @param[in] y_global a scalar or vector of classes.
     * @param[in] x design matrix
     * @param[in] beta weight vector
     * @param[in] cuts cutpoints vector
     * @param N_instances number of cases
     * @param N_attributes number of attributes
     * @param N_classes number of classes
     * @param is_y_vector 0 or 1 - whether y is a vector (alternatively it is a
     * scalar)
     * @param need_location_derivative interpreted as boolean - whether
     * location_derivative needs to be computed
     * @param need_cuts_derivative interpreted as boolean - whether
     * cuts_derivative needs to be computed
     */
    __kernel void ordered_logistic_glm(
        __global double* location_sum, __global double* logp_global,
        __global double* location_derivative, __global double* cuts_derivative,
        const __global int* y_global, const __global double* x,
        const __global double* beta, const __global double* cuts,
        const int N_instances, const int N_attributes, const int N_classes,
        const int is_y_vector, const int need_location_derivative,
        const int need_cuts_derivative) {
      const int gid = get_global_id(0);
      const int lid = get_local_id(0);
      const int lsize = get_local_size(0);
      const int wg_id = get_group_id(0);
      const int ngroups = get_num_groups(0);

      __local double local_storage[LOCAL_SIZE_];

      double logp = 0;
      double d1 = 0;
      double d2 = 0;
      double location = 0;
      int y;
      // Most calculations only happen for relevant data within next if.
      // Exceptions are reductions between threads that need barriers.
      if (gid < N_instances) {
        for (int i = 0, j = 0; i < N_attributes; i++, j += N_instances) {
          location += x[j + gid] * beta[i];
        }
        y = y_global[gid * is_y_vector];
        if (y < 1 || y > N_classes) {
          location = NAN;
        } else {
          const double cut_y1 = y == N_classes ? INFINITY : cuts[y - 1];
          const double cut_y2 = y == 1 ? -INFINITY : cuts[y - 2];
          const double cut1 = location - cut_y1;
          const double cut2 = location - cut_y2;

          if (y != N_classes) {
            logp -= log1p_exp(cut1);
          }
          if (y != 1) {
            logp -= log1p_exp(-cut2);
          }
          if (y != 1 && y != N_classes) {
            logp += log1m_exp(cut1 - cut2);
          }

          if (need_location_derivative || need_cuts_derivative) {
            double exp_cuts_diff = exp(cut_y2 - cut_y1);
            d1 = inv_logit(-cut2);
            d1 -= exp_cuts_diff / (exp_cuts_diff - 1);
            d2 = 1 / (1 - exp_cuts_diff);
            d2 -= inv_logit(-cut1);

            if (need_location_derivative) {
              location_derivative[gid] = d1 - d2;
            }
          }
        }
      }
      if (need_cuts_derivative) {
        for (int i = 0; i < N_classes - 1; i++) {
          local_storage[lid] = 0;
          if (gid < N_instances) {
            if (y - 1 == i) {
              local_storage[lid] = d2;
            } else if (y - 2 == i) {
              local_storage[lid] = -d1;
            }
          }
          // Sum cuts_derivative, calculated by different threads.
          // Since we can't sum between different work groups, we emit one
          // number per work group. These must be summed on CPU for final
          // result.
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
            cuts_derivative[(N_classes - 1) * wg_id + i] = local_storage[0];
          }
          barrier(CLK_LOCAL_MEM_FENCE);
        }
      }
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

      barrier(CLK_LOCAL_MEM_FENCE);
      local_storage[lid] = location;
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
        location_sum[wg_id] = local_storage[0];
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/ordered_logistic_glm_lpmf.hpp
 * ordered_logistic_glm() \endlink
 */
const kernel_cl<out_buffer, out_buffer, out_buffer, out_buffer, in_buffer,
                in_buffer, in_buffer, in_buffer, int, int, int, int, int, int>
    ordered_logistic_glm("ordered_logistic_glm",
                         {log1p_exp_device_function, log1m_exp_device_function,
                          inv_logit_device_function,
                          ordered_logistic_glm_kernel_code},
                         {{"REDUCTION_STEP_SIZE", 4}, {"LOCAL_SIZE_", 64}});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
