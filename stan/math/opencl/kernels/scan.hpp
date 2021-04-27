#ifndef STAN_MATH_OPENCL_KERNELS_SCAN_HPP
#define STAN_MATH_OPENCL_KERNELS_SCAN_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static const std::string scan_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     * Sums a batch of matrices. Buffer A contains
     * batch_size matrices of size rows x cols. All elements
     * at matching indices are summed up and stored to the
     * resulting matrix B.
     *
     * @param[out] B buffer of the result matrix.
     * @param[in] A buffer containing the entire batch.
     * @param rows Number of rows for a single matrix in the batch.
     * @param cols Number of cols for a single matrix in the batch.
     * @param batch_size Number of matrices in the batch.
     * @note Code is a <code>const char*</code> held in
     * <code>add_batch_kernel_code.</code>
     * This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void scan(__global SCAL *out_wgs, __global SCAL *out_threads,
                       __global SCAL *in, int size) {
      const int gid = get_global_id(0);
      const int lid = get_local_id(0);
      const int lsize = get_local_size(0);
      const int wg_id = get_group_id(0);
      const int gsize = get_global_size(0);

      int start = gid * size / gsize;
      int end = (gid + 1) * size / gsize;
      __local SCAL local_storage[LOCAL_SIZE_];

      SCAL acc = in[start];
      for (int i = start + 1; i < end; i++) {
        acc = OP(acc, in[i]);
      }
      out_threads[gid] = acc;
      local_storage[lid] = acc;
      barrier(CLK_LOCAL_MEM_FENCE);
      for (int step = lsize / REDUCTION_STEP_SIZE; step > 0;
           step /= REDUCTION_STEP_SIZE) {
        if (lid < step) {
          for (int i = 1; i < REDUCTION_STEP_SIZE; i++) {
            local_storage[lid]
                = OP(local_storage[lid], local_storage[lid + step * i]);
          }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      if (lid == 0) {
        out_wgs[wg_id] = local_storage[0];
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/add.hpp add_batch() \endlink
 */
const kernel_cl<out_buffer, in_buffer, int> scan_double_sum(
    "scan",
    {"#define SCAL double\n", "#define OP(x,y) ((x)+(y))\n", scan_kernel_code},
    {{"REDUCTION_STEP_SIZE", 4}, {"LOCAL_SIZE_", 64}});
}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
