#ifndef STAN_MATH_OPENCL_KERNELS_SCALAR_SUM_HPP
#define STAN_MATH_OPENCL_KERNELS_SCALAR_SUM_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static constexpr const char* scalar_sum_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     * Sums a buffer into the first element of another buffer using a single
     * work group. Each thread sums a strided slice of the input, then the
     * work group reduces the partial sums in local memory and one thread
     * writes the result.
     *
     * Must be launched with one work group whose size is a power of two no
     * larger than `LOCAL_SIZE_`.
     *
     * @param[in,out] out buffer whose first element receives the sum
     * @param[in] in values to sum
     * @param n number of values in `in`
     * @param offset value added to the sum
     * @param accumulate if nonzero the result is added to `out[0]`,
     * otherwise it overwrites `out[0]`
     */
    __kernel void scalar_sum(__global double* out, const __global double* in,
                             const int n, const double offset,
                             const int accumulate) {
      const int lid = get_local_id(0);
      const int lsize = get_local_size(0);
      __local double local_storage[LOCAL_SIZE_];
      double partial = 0;
      for (int i = lid; i < n; i += lsize) {
        partial += in[i];
      }
      local_storage[lid] = partial;
      barrier(CLK_LOCAL_MEM_FENCE);
      for (int step = lsize / 2; step > 0; step /= 2) {
        if (lid < step) {
          local_storage[lid] += local_storage[lid + step];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      if (lid == 0) {
        const double result = local_storage[0] + offset;
        if (accumulate) {
          out[0] += result;
        } else {
          out[0] = result;
        }
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/scalar_sum.hpp scalar_sum() \endlink
 */
const kernel_cl<in_out_buffer, in_buffer, int, double, int> scalar_sum(
    "scalar_sum", {scalar_sum_kernel_code}, {{"LOCAL_SIZE_", 256}});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
