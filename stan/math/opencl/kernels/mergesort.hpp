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
static const std::string mergesort_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     */
    void merge(__global SCAL* A, __global SCAL* B, __global SCAL* res,
               int A_size, int B_size) {
      int A_i = 0;
      int B_i = 0;
      if (B_i < B_size && A_i < A_size) {
        while (1) {
          if (A[A_i] < B[B_i]) {
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

    __kernel void merge_step(__global SCAL* output, __global SCAL* input,
                             int run_len, int size, int tasks) {
      int gid = get_global_id(0);
      int n_threads = get_global_size(0);
      for (int id = gid; id < tasks; id += n_threads) {
        int start_a = (2 * id) * run_len;
        int start_b = (2 * id + 1) * run_len;
        int end = min(size, (2 * id + 2) * run_len);
        //      printf("%d %d %d :%d %d %d \n", gid, n_threads, size, start_a,
        //      start_b,
        //             end);
        merge(input + start_a, input + start_b, output + start_a,
              start_b - start_a, end - start_b);
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/add.hpp add_batch() \endlink
 */
const kernel_cl<out_buffer, in_buffer, int, int, int> merge_step(
    "merge_step", {"#define SCAL double\n", mergesort_kernel_code});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
