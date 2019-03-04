#ifndef STAN_MATH_OPENCL_KERNELS_TRI_RECT_MULTIPLY_HPP
#define STAN_MATH_OPENCL_KERNELS_TRI_RECT_MULTIPLY_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static const char* tri_rect_multiply_kernel_code = STRINGIFY(
    // \endcond
    /**
     * OpenCL Matrix multiplication of a triangular and rectangular matrix 
     *
     * @param[in] A the lower/upper triangular matrix in matrix multiplication
     * @param[in] B the right matrix in matrix multiplication
     * @param[out] C the output matrix
     * @param[in] M Number of rows for matrix A
     * @param[in] N Number of rows for matrix B
     * @param[in] K Number of cols for matrix A and number of rows for matrix B
     * @param[in] lower_upper the triangularity of A (lower or upper)
     */
    __kernel void tri_rect_multiply(const __global double* A,
                                  const __global double* B, __global double* C,
                                  const int M, const int N, const int K,
                                  unsigned int lower_upper) {
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
        int do_mult = 0;
        // check if the of multiplication of the tiles must be performed
        // whole tile multiplies have to be performed even if some thread
        // are above/below the diagonal for lower/upper matrices. This is
        // due to the shared memory mechanism where all threads in the block
        // copy the data.
        // if A is lower tri and the smallest tile index in the thread
        // block is smaller than the row number, then perform the tile multiply
        do_mult |= lower_upper == LOWER &&
                   ((THREAD_BLOCK_SIZE * (tile_ind)) <= i);
        // if A is upper tri and the largest index in the thread block
        // is larger than the row number, then perform the tile multiply
        do_mult |= lower_upper == UPPER &&
                   ((THREAD_BLOCK_SIZE * (tile_ind+1)) >= i);
        if (do_mult) {
          const int tiled_i = THREAD_BLOCK_SIZE * tile_ind + thread_block_row;
          const int tiled_j = THREAD_BLOCK_SIZE * tile_ind + thread_block_col;
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
          // wait until all tile values are loaded to the local memory
          barrier(CLK_LOCAL_MEM_FENCE);
          for (int block_ind = 0; block_ind < THREAD_BLOCK_SIZE; block_ind++) {
            for (int w = 0; w < WORK_PER_THREAD; w++) {
                acc[w] += A_local[block_ind][thread_block_row]
                          * B_local[thread_block_col
                          * + w * THREAD_BLOCK_SIZE_COL][block_ind];
            }
          }
          barrier(CLK_LOCAL_MEM_FENCE);
        }
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
const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer,
                         int, int, int, TriangularViewCL>
    tri_rect_multiply("tri_rect_multiply", tri_rect_multiply_kernel_code,
                    {{"THREAD_BLOCK_SIZE", 32}, {"WORK_PER_THREAD", 8}});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
