#ifndef STAN_MATH_OPENCL_KERNELS_SCALAR_REDUCE_HPP
#define STAN_MATH_OPENCL_KERNELS_SCALAR_REDUCE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static constexpr const char* scalar_sum_op
    = "#define SCALAR_REDUCE_OP(a, b) ((a) + (b))\n"
      "#define SCALAR_REDUCE_INIT 0.0\n";
static constexpr const char* scalar_max_op
    = "#define SCALAR_REDUCE_OP(a, b) fmax((a), (b))\n"
      "#define SCALAR_REDUCE_INIT (-INFINITY)\n";
static constexpr const char* scalar_prod_op
    = "#define SCALAR_REDUCE_OP(a, b) ((a) * (b))\n"
      "#define SCALAR_REDUCE_INIT 1.0\n";

static constexpr const char* scalar_reduce_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     * Reduces a buffer into the first element of another buffer using a
     * single work group. Each thread reduces a strided slice of the input,
     * then the work group reduces the partial results in local memory and
     * one thread writes the result. The reduction is given by the macros
     * `SCALAR_REDUCE_OP(a, b)` and `SCALAR_REDUCE_INIT`.
     *
     * Must be launched with one work group whose size is a power of two no
     * larger than `LOCAL_SIZE_`.
     *
     * @param[in,out] out buffer whose first element receives the result
     * @param[in] in values to reduce
     * @param n number of values in `in`
     * @param offset value added to the result
     * @param accumulate if nonzero the result is combined with `out[0]`
     * using the reduction, otherwise it overwrites `out[0]`
     */
    __kernel void scalar_reduce(__global double* out, const __global double* in,
                                const int n, const double offset,
                                const int accumulate) {
      const int lid = get_local_id(0);
      const int lsize = get_local_size(0);
      __local double local_storage[LOCAL_SIZE_];
      double partial = SCALAR_REDUCE_INIT;
      for (int i = lid; i < n; i += lsize) {
        partial = SCALAR_REDUCE_OP(partial, in[i]);
      }
      local_storage[lid] = partial;
      barrier(CLK_LOCAL_MEM_FENCE);
      for (int step = lsize / 2; step > 0; step /= 2) {
        if (lid < step) {
          local_storage[lid]
              = SCALAR_REDUCE_OP(local_storage[lid], local_storage[lid + step]);
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      if (lid == 0) {
        const double result = local_storage[0] + offset;
        if (accumulate) {
          out[0] = SCALAR_REDUCE_OP(out[0], result);
        } else {
          out[0] = result;
        }
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * Sums a buffer into a device scalar. See the docs for \link
 * kernels/scalar_reduce.hpp scalar_reduce() \endlink
 */
const kernel_cl<in_out_buffer, in_buffer, int, double, int> scalar_sum(
    "scalar_reduce", {scalar_sum_op, scalar_reduce_kernel_code},
    {{"LOCAL_SIZE_", 256}});

/** \ingroup opencl_kernels
 * Finds the maximum of a buffer as a device scalar. See the docs for \link
 * kernels/scalar_reduce.hpp scalar_reduce() \endlink
 */
const kernel_cl<in_out_buffer, in_buffer, int, double, int> scalar_max(
    "scalar_reduce", {scalar_max_op, scalar_reduce_kernel_code},
    {{"LOCAL_SIZE_", 256}});

/** \ingroup opencl_kernels
 * Multiplies the elements of a buffer into a device scalar. See the docs for
 * \link kernels/scalar_reduce.hpp scalar_reduce() \endlink
 */
const kernel_cl<in_out_buffer, in_buffer, int, double, int> scalar_prod(
    "scalar_reduce", {scalar_prod_op, scalar_reduce_kernel_code},
    {{"LOCAL_SIZE_", 256}});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
