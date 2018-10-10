#ifndef STAN_MATH_GPU_KERNELS_LOWER_TRI_INVERSE_STEP3_HPP
#define STAN_MATH_GPU_KERNELS_LOWER_TRI_INVERSE_STEP3_HPP
#ifdef STAN_OPENCL

#include <stan/math/gpu/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
const char* lower_tri_inverse_step3_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Calculates the products to calculate parts of the lower_tri_inverse
     *
     * @param[in, out] A Input matrix that is being inverted.
     * @param[in] C Temporary matrix with the intermediate results.
     * @param A_rows Number of rows for A.
     * @param rows The number of rows in a single matrix of the batch
     * @param non_padded_rows Number of rows in A not used for padding.
     * @note Code is a <code>const char*</code> held in
     * <code>lower_tri_inverse_step3_kernel_code.</code>
     *  Used in math/gpu/lower_tri_inverse.hpp.
     *  This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void lower_tri_inverse_step3(
        __global read_write double* A, const __global double* temp,
        const int A_rows, const int rows, int non_padded_rows) {
      int result_matrix_id = get_global_id(2);
      int offset = result_matrix_id * rows * 2;
      // thread index inside the thread_block
      const int thread_block_row = get_local_id(0);
      const int thread_block_col = get_local_id(1);
      // global thread index
      const int i = THREAD_BLOCK_SIZE * get_group_id(0) + thread_block_row;
      const int j = THREAD_BLOCK_SIZE * get_group_id(1) + thread_block_col;

      // local memory
      __local double temp_local[THREAD_BLOCK_SIZE][THREAD_BLOCK_SIZE];
      __local double C1_local[THREAD_BLOCK_SIZE][THREAD_BLOCK_SIZE];

      double acc[WORK_PER_THREAD] = {0};

      const int num_tiles = (rows + THREAD_BLOCK_SIZE - 1) / THREAD_BLOCK_SIZE;
      // iterate over all tiles
      for (int tile_ind = 0; tile_ind < num_tiles; tile_ind++) {
        // each thread copies WORK_PER_THREAD values to the local
        // memory
        for (int w = 0; w < WORK_PER_THREAD; w++) {
          const int tiled_i = THREAD_BLOCK_SIZE * tile_ind + thread_block_row;
          const int tiled_j = THREAD_BLOCK_SIZE * tile_ind + thread_block_col;
          const int temp_global_col = tiled_j + w * THREAD_BLOCK_SIZE_COL;
          const int C1_global_col = offset + j + w * THREAD_BLOCK_SIZE_COL;
          const int C1_global_row = tiled_i + offset;
          if ((temp_global_col) < rows && (i) < rows) {
            temp_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                      [thread_block_row]
                = temp[result_matrix_id * rows * rows + temp_global_col * rows
                       + i];
          } else {
            temp_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                      [thread_block_row]
                = 0.0;
          }
          if (C1_global_col <= C1_global_row) {
            C1_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                    [thread_block_row]
                = A[C1_global_col * A_rows + C1_global_row];
          } else {
            C1_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                    [thread_block_row]
                = 0;
          }
        }
        // wait until all tile values are loaded to the local memory
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int block_ind = 0; block_ind < THREAD_BLOCK_SIZE; block_ind++) {
          for (int w = 0; w < WORK_PER_THREAD; w++) {
            acc[w] += temp_local[block_ind][thread_block_row]
                      * C1_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                                [block_ind];
          }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      // save the values
      for (int w = 0; w < WORK_PER_THREAD; w++) {
        // each thread saves WORK_PER_THREAD values
        A[(offset + j + w * THREAD_BLOCK_SIZE_COL) * A_rows + i + rows + offset]
            = -acc[w];
      }
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/matrix_multiply.hpp add() \endlink
 */
const local_range_kernel<cl::Buffer, cl::Buffer, int, int, int>
    lower_tri_inverse_step3("lower_tri_inverse_step3",
                            lower_tri_inverse_step3_kernel_code,
                            {{"THREAD_BLOCK_SIZE", 32},
                             {"WORK_PER_THREAD", 8}});
}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
