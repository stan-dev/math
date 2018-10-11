#ifndef STAN_MATH_GPU_KERNELS_LOWER_TRI_INVERSE_HPP
#define STAN_MATH_GPU_KERNELS_LOWER_TRI_INVERSE_HPP
#ifdef STAN_OPENCL

#include <stan/math/gpu/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
const char* lower_tri_inverse_step2_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Calculates the intermediate products in the lower_tri_inverse
     *
     * @param[in] A input matrix that is being inverted.
     * @param[out] temp output matrix with results of the batched matrix
     * multiplications
     * @param A_rows The number of rows for A.
     * @param rows The number of rows in a single matrix of the batch
     * @note Code is a <code>const char*</code> held in
     * <code>lower_tri_inverse_step2_kernel_code.</code>
     *  Used in math/gpu/lower_tri_inverse.hpp.
     *  This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void lower_tri_inverse_step2(__global double* A,
                                          __global double* temp,
                                          const int A_rows, const int rows) {
      int result_matrix_id = get_global_id(2);
      int offset = result_matrix_id * rows * 2;
      // thread index inside the thread_block
      const int thread_block_row = get_local_id(0);
      const int thread_block_col = get_local_id(1);
      // global thread index
      const int global_thread_row
          = THREAD_BLOCK_SIZE * get_group_id(0) + thread_block_row;
      const int global_thread_col
          = THREAD_BLOCK_SIZE * get_group_id(1) + thread_block_col;

      // local memory
      __local double C2_local[THREAD_BLOCK_SIZE][THREAD_BLOCK_SIZE];
      __local double A3_local[THREAD_BLOCK_SIZE][THREAD_BLOCK_SIZE];

      double acc[WORK_PER_THREAD] = {0};

      const int num_tiles = (rows + THREAD_BLOCK_SIZE - 1) / THREAD_BLOCK_SIZE;
      // iterate over all tiles
      for (int tile_ind = 0; tile_ind < num_tiles; tile_ind++) {
        // each thread copies WORK_PER_THREAD values to the local
        // memory
        for (int w = 0; w < WORK_PER_THREAD; w++) {
          const int tiled_i = THREAD_BLOCK_SIZE * tile_ind + thread_block_row;
          const int tiled_j = THREAD_BLOCK_SIZE * tile_ind + thread_block_col;

          const int C2_global_col
              = offset + rows + tiled_j + w * THREAD_BLOCK_SIZE_COL;
          const int C2_global_row = offset + global_thread_row + rows;
          const int A3_global_col
              = offset + global_thread_col + w * THREAD_BLOCK_SIZE_COL;
          const int A3_global_row = tiled_i + rows + offset;
          const int local_col = thread_block_col + w * THREAD_BLOCK_SIZE_COL;
          const int local_row = thread_block_row;
          if (C2_global_col <= C2_global_row) {
            C2_local[local_col][local_row]
                = A[C2_global_col * A_rows + C2_global_row];
          } else {
            C2_local[local_col][local_row] = 0;
          }
          A3_local[local_col][local_row]
              = A[A3_global_col * A_rows + A3_global_row];
        }
        // wait until all tile values are loaded to the local memory
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int block_ind = 0; block_ind < THREAD_BLOCK_SIZE; block_ind++) {
          for (int w = 0; w < WORK_PER_THREAD; w++) {
            const int local_col = thread_block_col + w * THREAD_BLOCK_SIZE_COL;
            const int local_row = thread_block_row;
            acc[w] += C2_local[block_ind][local_row]
                      * A3_local[local_col][block_ind];
          }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      const int batch_offset = result_matrix_id * rows * rows;
      const int temp_global_row = global_thread_row;
      // save the values
      for (int w = 0; w < WORK_PER_THREAD; w++) {
        // each thread saves WORK_PER_THREAD values
        const int temp_global_col
            = global_thread_col + w * THREAD_BLOCK_SIZE_COL;

        temp[batch_offset + temp_global_col * rows + temp_global_row] = acc[w];
      }
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/matrix_multiply.hpp add() \endlink
 */
const local_range_kernel<cl::Buffer, cl::Buffer, int, int>
    lower_tri_inverse_step2("lower_tri_inverse_step2",
                            lower_tri_inverse_step2_kernel_code,
                            {{"THREAD_BLOCK_SIZE", 32},
                             {"WORK_PER_THREAD", 8}});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
