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

      int start = (int)((long)gid * size / gsize);
      int end = (int)((long)(gid + 1) * size / gsize);
      __local SCAL local_storage[LOCAL_SIZE_];

      SCAL acc = 0;
      if (start != end) {
        acc = in[start];
        for (int i = start + 1; i < end; i++) {
          acc += in[i];
        }
      }
      for (int step = 1; step < lsize; step *= REDUCTION_STEP_SIZE) {
        local_storage[lid] = acc;
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int i = 1; i < REDUCTION_STEP_SIZE && step * i <= lid; i++) {
          acc += local_storage[lid - step * i];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      out_threads[gid] = acc;
      if (lid == LOCAL_SIZE_ - 1) {
        out_wgs[wg_id] = acc;
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/add.hpp add_batch() \endlink
 */
const kernel_cl<out_buffer, out_buffer, in_buffer, int> scan_double_sum(
    "scan", {"#define SCAL double\n", scan_kernel_code},
    {{"REDUCTION_STEP_SIZE", 4}, {"LOCAL_SIZE_", 16}});

// \cond
static const std::string scan2_kernel_code = STRINGIFY(
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
    __kernel void scan2(__global SCAL *data, int size) {
      const int gid = get_global_id(0);
      const int gsize = get_global_size(0);

      int start = (int)((long)gid * size / gsize);
      int end = (int)((long)(gid + 1) * size / gsize);
      __local SCAL local_storage[LOCAL_SIZE_];

      SCAL acc;
      if (start == end) {
        acc = 0;
      } else {
        acc = data[start];
        for (int i = start + 1; i < end; i++) {
          acc += data[i];
        }
      }
      local_storage[gid] = acc;
      barrier(CLK_LOCAL_MEM_FENCE);
      for (int step = 1; step < gsize; step *= REDUCTION_STEP_SIZE) {
        for (int i = 1; i < REDUCTION_STEP_SIZE && step * i <= gid; i++) {
          acc += local_storage[gid - step * i];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        local_storage[gid] = acc;
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      if (start != end) {
        if (gid == 0) {
          acc = 0;
        } else {
          acc = local_storage[gid - 1];
        }
        for (int i = start; i < end; i++) {
          acc += data[i];
          data[i] = acc;
        }
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/add.hpp add_batch() \endlink
 */
const kernel_cl<in_out_buffer, int> scan2_double_sum(
    "scan2", {"#define SCAL double\n", scan2_kernel_code},
    {{"REDUCTION_STEP_SIZE", 4}, {"LOCAL_SIZE_", 1024}});

// \cond
static const std::string scan3_kernel_code = STRINGIFY(
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
    __kernel void scan3(__global SCAL *out, __global SCAL *in_data,
                        __global SCAL *in_threads, __global SCAL *in_wgs,
                        int size) {
      const int gid = get_global_id(0);
      const int lid = get_local_id(0);
      const int lsize = get_local_size(0);
      const int wg_id = get_group_id(0);
      const int gsize = get_global_size(0);

      int start = (int)((long)gid * size / gsize);
      int end = (int)((long)(gid + 1) * size / gsize);
      __local SCAL local_storage[LOCAL_SIZE_];

      SCAL acc = 0;
      if (wg_id != 0) {
        acc = in_wgs[wg_id - 1];
      }
      if (lid != 0) {
        acc += in_threads[gid - 1];
      }
      for (int i = start; i < end; i++) {
        acc += in_data[i];
        out[i] = acc;
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/add.hpp add_batch() \endlink
 */
const kernel_cl<out_buffer, in_buffer, in_buffer, in_buffer, int>
    scan3_double_sum("scan3", {"#define SCAL double\n", scan3_kernel_code},
                     {{"REDUCTION_STEP_SIZE", 4}, {"LOCAL_SIZE_", 16}});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
