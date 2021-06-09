#ifndef STAN_MATH_OPENCL_KERNELS_CUMULATIVE_SUM_HPP
#define STAN_MATH_OPENCL_KERNELS_CUMULATIVE_SUM_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static const char *cumulative_sum1_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     * First kernel of the cumulative sum implementation. Each thread sums the
     * assigned elements and threads within same work group add their results
     * together.
     *
     * @param[out] out_wgs results from each work group
     * @param[out] out_threads results for each thread
     * @param[in] in input data
     * @param size size number of elements in the input
     */
    __kernel void cumulative_sum1(__global SCAL *out_wgs,
                                  __global SCAL *out_threads, __global SCAL *in,
                                  int size) {
      const int gid = get_global_id(0);
      const int lid = get_local_id(0);
      const int lsize = get_local_size(0);
      const int wg_id = get_group_id(0);
      const int gsize = get_global_size(0);

      int start = (int)((long)gid * size / gsize);      // NOLINT
      int end = (int)((long)(gid + 1) * size / gsize);  // NOLINT
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

// \cond
static const char *cumulative_sum2_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     * Second kernel of the cumulative sum implementation. Calculates prefix sum
     * of given data in place using a single work group (must be run with a
     * single work group).
     *
     * @param[in, out] data data to calculate cumulative sum of
     * @param size size number of elements in the input
     */
    __kernel void cumulative_sum2(__global SCAL *data, int size) {
      const int gid = get_global_id(0);
      const int gsize = get_global_size(0);

      int start = (int)((long)gid * size / gsize);      // NOLINT
      int end = (int)((long)(gid + 1) * size / gsize);  // NOLINT
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

// \cond
static const char *cumulative_sum3_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     * Third kernel of the cumulative sum implementation. Given sums of threads
     * and cumulative sum of those calculates cumulative sum of given array.
     * Must be run with the same number of threads and work groups as the first
     * cumulative sum kernel.
     *
     * @param[out] out cumulatively summed input
     * @param[out] in_data input data
     * @param[in] in_threads summed results from each thread from the first
     * kernel
     * @param[in] in_wgs cumulatively summed results from each work group
     * (calculated by previous two kernels)
     * @param size size number of elements in the input
     */
    __kernel void cumulative_sum3(__global SCAL *out, __global SCAL *in_data,
                                  __global SCAL *in_threads,
                                  __global SCAL *in_wgs, int size) {
      const int gid = get_global_id(0);
      const int lid = get_local_id(0);
      const int lsize = get_local_size(0);
      const int wg_id = get_group_id(0);
      const int gsize = get_global_size(0);

      int start = (int)((long)gid * size / gsize);      // NOLINT
      int end = (int)((long)(gid + 1) * size / gsize);  // NOLINT
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

/**
 * struct containing cumulative_sum kernels, grouped by scalar type.
 */
template <typename Scalar, typename = void>
struct cumulative_sum {};

template <typename T>
struct cumulative_sum<double, T> {
  static const kernel_cl<out_buffer, out_buffer, in_buffer, int> kernel1;
  static const kernel_cl<in_out_buffer, int> kernel2;
  static const kernel_cl<out_buffer, in_buffer, in_buffer, in_buffer, int>
      kernel3;
};
template <typename T>
struct cumulative_sum<int, T> {
  static const kernel_cl<out_buffer, out_buffer, in_buffer, int> kernel1;
  static const kernel_cl<in_out_buffer, int> kernel2;
  static const kernel_cl<out_buffer, in_buffer, in_buffer, in_buffer, int>
      kernel3;
};

template <typename T>
const kernel_cl<out_buffer, out_buffer, in_buffer, int>
    cumulative_sum<double, T>::kernel1("cumulative_sum1",
                                       {"#define SCAL double\n",
                                        cumulative_sum1_kernel_code},
                                       {{"REDUCTION_STEP_SIZE", 4},
                                        {"LOCAL_SIZE_", 16}});
template <typename T>
const kernel_cl<out_buffer, out_buffer, in_buffer, int>
    cumulative_sum<int, T>::kernel1(
        "cumulative_sum1", {"#define SCAL int\n", cumulative_sum1_kernel_code},
        {{"REDUCTION_STEP_SIZE", 4}, {"LOCAL_SIZE_", 16}});

template <typename T>
const kernel_cl<in_out_buffer, int> cumulative_sum<double, T>::kernel2(
    "cumulative_sum2", {"#define SCAL double\n", cumulative_sum2_kernel_code},
    {{"REDUCTION_STEP_SIZE", 4}, {"LOCAL_SIZE_", 1024}});
template <typename T>
const kernel_cl<in_out_buffer, int> cumulative_sum<int, T>::kernel2(
    "cumulative_sum2", {"#define SCAL int\n", cumulative_sum2_kernel_code},
    {{"REDUCTION_STEP_SIZE", 4}, {"LOCAL_SIZE_", 1024}});

template <typename T>
const kernel_cl<out_buffer, in_buffer, in_buffer, in_buffer, int>
    cumulative_sum<double, T>::kernel3("cumulative_sum3",
                                       {"#define SCAL double\n",
                                        cumulative_sum3_kernel_code},
                                       {{"REDUCTION_STEP_SIZE", 4},
                                        {"LOCAL_SIZE_", 16}});
template <typename T>
const kernel_cl<out_buffer, in_buffer, in_buffer, in_buffer, int>
    cumulative_sum<int, T>::kernel3(
        "cumulative_sum3", {"#define SCAL int\n", cumulative_sum3_kernel_code},
        {{"REDUCTION_STEP_SIZE", 4}, {"LOCAL_SIZE_", 16}});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
