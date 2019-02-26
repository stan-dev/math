#ifndef STAN_MATH_GPU_KERNELS_LOWER_TRI_RECT_MULTIPLY_HPP
#define STAN_MATH_GPU_KERNELS_LOWER_TRI_RECT_MULTIPLY_HPP
#ifdef STAN_OPENCL

#include <stan/math/gpu/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static const char* lower_tri_rect_multiply_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Matrix multiplication of a lower triangular and rectengular matrix 
     * on the GPU
     *
     * @param[in] A the left lower triangular matrix in matrix multiplication
     * @param[in] B the right matrix in matrix multiplication
     * @param[out] C the output matrix
     * @param[in] M Number of rows for matrix A
     * @param[in] N Number of rows for matrix B
     * @param[in] K Number of cols for matrix A and number of rows for matrix B
     */
    __kernel void lower_tri_rect_multiply(const __global double* A,
                                  const __global double* B, __global double* C,
                                  const int M, const int N, const int K) {
      // thread index inside the thread_block
      const int thread_block_row = get_local_id(0);
      const int thread_block_col = get_local_id(1);
      // global thread index
      const int i = THREAD_BLOCK_SIZE * get_group_id(0) + thread_block_row;
      const int j = THREAD_BLOCK_SIZE * get_group_id(1) + thread_block_col;

      // indexes that determine the last indexes that need to compute
      // in order to remove the unnecesary multiplications in the special
      // multiplication of A*A^T
      const int jMin = THREAD_BLOCK_SIZE * get_group_id(1);
      const int iMax = THREAD_BLOCK_SIZE * get_group_id(0) + get_local_size(0);
      // local memory
      __local double A_local[THREAD_BLOCK_SIZE][THREAD_BLOCK_SIZE];
      __local double B_local[THREAD_BLOCK_SIZE][THREAD_BLOCK_SIZE];

      double acc[WORK_PER_THREAD];
      for (int w = 0; w < WORK_PER_THREAD; w++) {
        acc[w] = 0.0;
      }

      const int num_tiles = (K + THREAD_BLOCK_SIZE - 1) / THREAD_BLOCK_SIZE;
      // iterate over all tiles
      for (int tile_ind = 0; tile_ind < num_tiles; tile_ind++) {
        const int tiled_i = THREAD_BLOCK_SIZE * tile_ind + thread_block_row;
        const int tiled_j = THREAD_BLOCK_SIZE * tile_ind + thread_block_col;
        // if the data needs to be loaded to local memory
        if (jMin >= iMax) {
          // each thread copies WORK_PER_THREAD values to the local
          // memory
          for (int w = 0; w < WORK_PER_THREAD; w++) {
            A_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                  [thread_block_row]
                = A[(tiled_j + w * THREAD_BLOCK_SIZE_COL) * M + i];

            B_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                  [thread_block_row]
                = B[(j + w * THREAD_BLOCK_SIZE_COL) * K + tiled_i];
          }
        }
        // wait until all tile values are loaded to the local memory
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int block_ind = 0; block_ind < THREAD_BLOCK_SIZE; block_ind++) {
          for (int w = 0; w < WORK_PER_THREAD; w++) {
            if (jMin >= iMax) {
            if ((THREAD_BLOCK_SIZE * tile_ind + block_ind) <= i) {
              acc[w] += A_local[block_ind][thread_block_row]
                        * B_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                                [block_ind];
            }
            }
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
 * See the docs for \link kernels/lower_tri_rect_multiply.hpp add() \endlink
 */
const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, int, int, int>
    lower_tri_rect_multiply("lower_tri_rect_multiply", lower_tri_rect_multiply_kernel_code,
                    {{"THREAD_BLOCK_SIZE", 32}, {"WORK_PER_THREAD", 8}});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
