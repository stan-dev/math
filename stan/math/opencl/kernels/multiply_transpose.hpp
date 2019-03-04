#ifndef STAN_MATH_OPENCL_KERNELS_MULTIPLY_TRANSPOSE_HPP
#define STAN_MATH_OPENCL_KERNELS_MULTIPLY_TRANSPOSE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static const char* multiply_transpose_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Matrix multiplication of the form A*A^T on the OpenCL device
     *
     * @param[in] A matrix A
     * @param[out] B the output matrix
     * @param[in] M Number of rows for matrix A
     * @param[in] N Number of cols for matrix A and the number of rows for
     * matrix A^T
     */
    __kernel void multiply_transpose(const __global double* A,
                                     __global double* B, const int M,
                                     const int N) {
      // thread index inside the thread block
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

      const int numTiles = (N + THREAD_BLOCK_SIZE - 1) / THREAD_BLOCK_SIZE;
      // iterate over all tiles
      for (int tile_ind = 0; tile_ind < numTiles; tile_ind++) {
        // in each tile
        const int tiled_i = THREAD_BLOCK_SIZE * tile_ind + thread_block_row;
        const int tiled_j = THREAD_BLOCK_SIZE * tile_ind + thread_block_col;
        // if the data needs to be loaded to local memory
        if (jMin <= iMax) {
          // each thread copies WORK_PER_THREAD values to the
          // local memory
          for (int w = 0; w < WORK_PER_THREAD; w++) {
            A_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                   [thread_block_row]
                = A[i + (tiled_j + w * THREAD_BLOCK_SIZE_COL) * M];
            B_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                   [thread_block_row]
                = A[(j + w * THREAD_BLOCK_SIZE_COL) + tiled_i * M];
          }
        }
        // wait till all tile values are loaded to the local memory
        barrier(CLK_LOCAL_MEM_FENCE);
        // multiply the tile products
        for (int block_ind = 0; block_ind < THREAD_BLOCK_SIZE; block_ind++) {
          // each thread multiplies WORK_PER_THREAD values
          for (int w = 0; w < WORK_PER_THREAD; w++) {
            if (jMin <= iMax) {
              if ((j + w * THREAD_BLOCK_SIZE_COL) <= i) {
                acc[w] += A_local[block_ind][thread_block_row]
                          * B_local[thread_block_col
                                    + w * THREAD_BLOCK_SIZE_COL][block_ind];
              }
            }
          }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      // save the values
      for (int w = 0; w < WORK_PER_THREAD; w++) {
        // each thread saves WORK_PER_THREAD values
        if ((j + w * THREAD_BLOCK_SIZE_COL) <= i) {
          B[i + (j + w * THREAD_BLOCK_SIZE_COL) * M] = acc[w];
          B[(j + w * THREAD_BLOCK_SIZE_COL) + i * M] = acc[w];
        }
      }
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/multiply_transpose.hpp add() \endlink
 */
const local_range_kernel<cl::Buffer, cl::Buffer, int, int> multiply_transpose(
    "multiply_transpose", multiply_transpose_kernel_code,
    {{"THREAD_BLOCK_SIZE", 32}, {"WORK_PER_THREAD", 4}});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
