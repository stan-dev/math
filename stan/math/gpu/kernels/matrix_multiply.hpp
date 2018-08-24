#ifndef STAN_MATH_GPU_KERNELS_MATRIX_MULTIPLY_HPP
#define STAN_MATH_GPU_KERNELS_MATRIX_MULTIPLY_HPP
#ifdef STAN_OPENCL

#include <stan/math/gpu/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
const char* matrix_multiply_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Matrix multiplication on the GPU
     *
     * @param[in] A the left matrix in matrix multiplication
     * @param[in] B the right matrix in matrix multiplication
     * @param[out] C the output matrix
     * @param[in] M Number of rows for matrix A
     * @param[in] N Number of rows for matrix B
     * @param[in] K Number of cols for matrix A and number of rows for matrix B
     */
    __kernel void matrix_multiply(const __global double* A,
                                  const __global double* B, __global double* C,
                                  const int M, const int N, const int K) {
      // thread index inside the workgroup
      const int workgroup_row = get_local_id(0);
      const int workgroup_col = get_local_id(1);
      // global thread index
      const int i = WORKGROUP_SIZE_MULT * get_group_id(0) + workgroup_row;
      const int j = WORKGROUP_SIZE_MULT * get_group_id(1) + workgroup_col;

      // local memory
      __local double A_local[WORKGROUP_SIZE_MULT][WORKGROUP_SIZE_MULT];
      __local double B_local[WORKGROUP_SIZE_MULT][WORKGROUP_SIZE_MULT];

      double acc[WORK_PER_THREAD_MULT];
      for (int w = 0; w < WORK_PER_THREAD_MULT; w++) {
        acc[w] = 0.0;
      }

      const int num_tiles = (K + WORKGROUP_SIZE_MULT - 1) / WORKGROUP_SIZE_MULT;
      // iterate over all tiles
      for (int t = 0; t < num_tiles; t++) {
        // each thread copies WORK_PER_THREAD_MULT values to the local
        // memory
        for (int w = 0; w < WORK_PER_THREAD_MULT; w++) {
          const int tiled_i = WORKGROUP_SIZE_MULT * t + workgroup_row;
          const int tiled_j = WORKGROUP_SIZE_MULT * t + workgroup_col;

          A_local[workgroup_col + w * WORKGROUP_SIZE_MULT_COL][workgroup_row]
              = A[(tiled_j + w * WORKGROUP_SIZE_MULT_COL) * M + i];

          B_local[workgroup_col + w * WORKGROUP_SIZE_MULT_COL][workgroup_row]
              = B[(j + w * WORKGROUP_SIZE_MULT_COL) * K + tiled_i];
        }
        // wait until all tile values are loaded to the local memory
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int k = 0; k < WORKGROUP_SIZE_MULT; k++) {
          for (int w = 0; w < WORK_PER_THREAD_MULT; w++) {
            acc[w] += A_local[k][workgroup_row]
                      * B_local[workgroup_col + w * WORKGROUP_SIZE_MULT_COL][k];
          }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      // save the values
      for (int w = 0; w < WORK_PER_THREAD_MULT; w++) {
        // each thread saves WORK_PER_THREAD_MULT_SELF_TRANS values
        C[(j + w * WORKGROUP_SIZE_MULT_COL) * M + i] = acc[w];
      }
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/matrix_multiply.hpp add() \endlink
 */
const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, int, int, int>
    matrix_multiply("matrix_multiply", matrix_multiply_kernel_code);

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
