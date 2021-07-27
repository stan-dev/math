#ifndef STAN_MATH_OPENCL_KERNELS_MERGESORT_HPP
#define STAN_MATH_OPENCL_KERNELS_MERGESORT_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static const char* mergesort_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Merges two sorted runs into a single sorted run of combined length.
     * @param A first sorted run
     * @param B second sorted run
     * @param[out] res result
     * @param A_size number of elements in `A`
     * @param B_size number of elements in `B`
     */
    void merge(__global SCAL* A, __global SCAL* B, __global SCAL* res,
               int A_size, int B_size) {
      int A_i = 0;
      int B_i = 0;
      if (B_i < B_size && A_i < A_size) {
        while (1) {
          if (A[A_i] CMP B[B_i]) {
            res[A_i + B_i] = A[A_i];
            A_i++;
            if (A_i >= A_size) {
              break;
            }
          } else {
            res[A_i + B_i] = B[B_i];
            B_i++;
            if (B_i >= B_size) {
              break;
            }
          }
        }
      }
      for (; A_i < A_size; A_i++) {
        res[A_i + B_i] = A[A_i];
      }
      for (; B_i < B_size; B_i++) {
        res[A_i + B_i] = B[B_i];
      }
    }

    /**
     * Searches for the index of the element that is larger than or equal to
     * given `value` in given range.
     * @param input data to search in
     * @param start start of the sorted range to search
     * @param end end of the sorted range to search
     * @param value value to search for
     */
    int binary_search(__global SCAL* input, int start, int end, SCAL value) {
      while (end > start) {
        int mid = start + (end - start) / 2;
        if (value CMP input[mid]) {
          end = mid;
        } else {
          start = mid + 1;
        }
      }
      return start;
    }

    /** \ingroup opencl_kernels
     *  Merges sorted runs into longer sorted runs.
     *  @param[out] output result
     *  @param input input data
     *  @param run_len lengths of sorted runs in the `input`
     *  @param size size of the input and output
     *  @param tasks number of merges that need to be executed - same as the
     * number of sorted runs in the output
     */
    __kernel void merge_step(__global SCAL* output, __global SCAL* input,
                             int run_len, int size, int tasks) {
      int gid = get_global_id(0);
      int n_threads = get_global_size(0);
      if (tasks >= n_threads) {
        // divide tasks between threads
        for (int task_id = gid; task_id < tasks; task_id += n_threads) {
          int start_a = (2 * task_id) * run_len;
          int start_b = min(size, (2 * task_id + 1) * run_len);
          int end = min(size, (2 * task_id + 2) * run_len);
          merge(input + start_a, input + start_b, output + start_a,
                start_b - start_a, end - start_b);
        }
      } else {
        // divide threads between tasks
        int task_id = gid % tasks;
        int threads_with_task_id = n_threads / tasks;
        threads_with_task_id
            += n_threads - threads_with_task_id * tasks > task_id;
        int id_in_task = gid / tasks;

        int start_a = (2 * task_id) * run_len;
        int start_b = min(size, (2 * task_id + 1) * run_len);
        int end = min(size, (2 * task_id + 2) * run_len);

        // divide a task between threads working on it
        int my_start_a = start_a
                         + (start_b - start_a) * (long)id_in_task  // NOLINT
                               / threads_with_task_id;
        int my_end_a = start_a
                       + (start_b - start_a) * (long)(id_in_task + 1)  // NOLINT
                             / threads_with_task_id;
        int my_start_b = id_in_task == 0 ? start_b
                                         : binary_search(input, start_b, end,
                                                         input[my_start_a]);
        int my_end_b
            = id_in_task == threads_with_task_id - 1
                  ? end
                  : binary_search(input, start_b, end, input[my_end_a]);

        int output_start
            = start_a + my_start_a - start_a + my_start_b - start_b;

        merge(input + my_start_a, input + my_start_b, output + output_start,
              my_end_a - my_start_a, my_end_b - my_start_b);
      }
    }
    // \cond
);
// \endcond

/**
 * struct containing sort_asc kernels, grouped by scalar type.
 */
template <typename Scalar, typename = void>
struct sort_asc {};

template <typename T>
struct sort_asc<double, T> {
  static const kernel_cl<out_buffer, in_buffer, int, int, int> merge_step;
};
template <typename T>
const kernel_cl<out_buffer, in_buffer, int, int, int>
    sort_asc<double, T>::merge_step("merge_step",
                                    {"#define SCAL double\n", "#define CMP <\n",
                                     mergesort_kernel_code});

template <typename T>
struct sort_asc<int, T> {
  static const kernel_cl<out_buffer, in_buffer, int, int, int> merge_step;
};
template <typename T>
const kernel_cl<out_buffer, in_buffer, int, int, int>
    sort_asc<int, T>::merge_step("merge_step",
                                 {"#define SCAL int\n", "#define CMP <\n",
                                  mergesort_kernel_code});

/**
 * struct containing sort_desc kernels, grouped by scalar type.
 */
template <typename Scalar, typename = void>
struct sort_desc {};

template <typename T>
struct sort_desc<double, T> {
  static const kernel_cl<out_buffer, in_buffer, int, int, int> merge_step;
};
template <typename T>
const kernel_cl<out_buffer, in_buffer, int, int, int>
    sort_desc<double, T>::merge_step("merge_step", {"#define SCAL double\n",
                                                    "#define CMP >\n",
                                                    mergesort_kernel_code});

template <typename T>
struct sort_desc<int, T> {
  static const kernel_cl<out_buffer, in_buffer, int, int, int> merge_step;
};
template <typename T>
const kernel_cl<out_buffer, in_buffer, int, int, int>
    sort_desc<int, T>::merge_step("merge_step",
                                  {"#define SCAL int\n", "#define CMP >\n",
                                   mergesort_kernel_code});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
