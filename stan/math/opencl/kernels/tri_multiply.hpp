#ifndef STAN_MATH_OPENCL_KERNELS_TRI_MULTIPLY_HPP
#define STAN_MATH_OPENCL_KERNELS_TRI_MULTIPLY_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static const char* tri_multiply_kernel_code = STRINGIFY(
    // \endcond
    /**
     * OpenCL Matrix multiplication of a triangular and rectangular matrix
     *
     * @param[in] A the lower/upper triangular matrix in matrix multiplication
     * @param[in] B the right matrix in matrix multiplication
     * @param[out] C the output matrix
     * @param[in] M Number of rows for matrix A
     * @param[in] N Number of cols for matrix B
     * @param[in] K Number of cols for matrix A and number of rows for matrix B
     * @param[in] lower_upper_A the triangularity of A (lower, upper or none)
     * @param[in] lower_upper_B the triangularity of B (lower, upper or none)
     */
    __kernel void tri_multiply(
        const __global double* A, const __global double* B, __global double* C,
        const int M, const int N, const int K, unsigned int lower_upper_A,
        unsigned int lower_upper_B) {
      // thread index inside the thread_block
      const int thread_block_row = get_local_id(0);
      const int thread_block_col = get_local_id(1);
      // global thread index
      const int i = THREAD_BLOCK_SIZE * get_group_id(0) + thread_block_row;
      const int j = THREAD_BLOCK_SIZE * get_group_id(1) + thread_block_col;

      // local memory
      __local double A_local[THREAD_BLOCK_SIZE][THREAD_BLOCK_SIZE];
      __local double B_local[THREAD_BLOCK_SIZE][THREAD_BLOCK_SIZE];

      double acc[WORK_PER_THREAD];
      for (int w = 0; w < WORK_PER_THREAD; w++) {
        acc[w] = 0.0;
      }

      const int num_tiles = (K + THREAD_BLOCK_SIZE - 1) / THREAD_BLOCK_SIZE;

      const int is_A_lower = lower_upper_A == LOWER;
      const int is_A_upper = lower_upper_A == UPPER;
      const int is_B_lower = lower_upper_B == LOWER;
      const int is_B_upper = lower_upper_B == UPPER;
      // This kernel is based on the kernel for the regular
      // matrix multiply where the multiplitcation is split
      // in tiles.
      // In cases where a matrix is lower triangular
      // its not necessary to multiply the elements
      // over the diagonal, therefore those tiles
      // in the matrix multiply can be skipped.
      // With upper triangular matrices we dont need
      // to multiply the elements under the diagonal,
      // so those tiles can be skipped.
      // The following code determines the start and
      // end tile based on triangularity of the input matrices
      // If no matrices are triangular the starting tile
      // is 0 and the end tile is num_tiles-1 which
      // is the same as in \link kernels/multiply.hpp multiply()
      const int end_tile_A
          = is_A_lower ? (i / THREAD_BLOCK_SIZE) : (num_tiles - 1);
      const int end_tile_B
          = is_B_upper ? (j / THREAD_BLOCK_SIZE) : (num_tiles - 1);
      const int start_tile_A = is_A_upper ? (i / THREAD_BLOCK_SIZE) : 0;
      const int start_tile_B = is_B_lower ? (j / THREAD_BLOCK_SIZE) : 0;
      const int start_tile = max(start_tile_A, start_tile_B);
      const int end_tile = min(end_tile_A, end_tile_B);

      for (int tile_idx = start_tile; tile_idx <= end_tile; tile_idx++) {
        const int tiled_i = THREAD_BLOCK_SIZE * tile_idx + thread_block_row;
        const int tiled_j = THREAD_BLOCK_SIZE * tile_idx + thread_block_col;
        // each thread copies WORK_PER_THREAD values to the local
        // memory
        for (int w = 0; w < WORK_PER_THREAD; w++) {
          // For the tiles on the diagonal we can ignore the values over
          // the diagonal if the matrix is lower triangular or under
          // the diagonal if the matrix is upper triangular
          if ((is_A_lower && ((tiled_j + w * THREAD_BLOCK_SIZE_COL) <= i))
              || (is_A_upper && ((tiled_j + w * THREAD_BLOCK_SIZE_COL) >= i))
              || (!is_A_lower && !is_A_upper)) {
            A_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                   [thread_block_row]
                = A[(tiled_j + w * THREAD_BLOCK_SIZE_COL) * M + i];
          } else {
            A_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                   [thread_block_row]
                = 0.0;
          }
          if ((is_B_lower && ((j + w * THREAD_BLOCK_SIZE_COL) <= tiled_i))
              || (is_B_upper && ((j + w * THREAD_BLOCK_SIZE_COL) >= tiled_i))
              || (!is_B_lower && !is_B_upper)) {
            B_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                   [thread_block_row]
                = B[(j + w * THREAD_BLOCK_SIZE_COL) * K + tiled_i];
          } else {
            B_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                   [thread_block_row]
                = 0.0;
          }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int block_idx = 0; block_idx < THREAD_BLOCK_SIZE; block_idx++) {
          for (int w = 0; w < WORK_PER_THREAD; w++) {
            acc[w] += A_local[block_idx][thread_block_row]
                      * B_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                               [block_idx];
          }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      // save the values
      for (int w = 0; w < WORK_PER_THREAD; w++) {
        // each thread saves WORK_PER_THREAD values
        C[(j + w * THREAD_BLOCK_SIZE_COL) * M + i] = acc[w];
      }
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/tri_multiply.hpp tri_multiply() \endlink
 */
const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, int, int, int,
                         TriangularViewCL, TriangularViewCL>
    tri_multiply("tri_multiply", tri_multiply_kernel_code,
                 {{"THREAD_BLOCK_SIZE", 32}, {"WORK_PER_THREAD", 8}});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
