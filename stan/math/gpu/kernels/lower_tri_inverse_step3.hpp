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
     * @param[in, out] A input matrix that is being inverted
     * @param[in] C tge temporary matrix with the intermediate results
     * @param M The number of rows for A.
     * @param K The number of elements to multiply
     * @param temp_rows The number of rows in C
     * @note Code is a <code>const char*</code> held in
     * <code>lower_tri_inverse_step3_kernel_code.</code>
     *  Used in math/gpu/lower_tri_inverse.hpp.
     *  This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void lower_tri_inverse_step3(
        __global read_write double* A, const __global read_only double* C,
        const read_only int M, const read_only int K,
        const read_only int temp_rows, int non_padded_rows) {
      int t = get_global_id(2);
      int offset = t * temp_rows * 2;
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
      // iterate over all tiles
      for (int tile_ind = 0; tile_ind < num_tiles; tile_ind++) {
        // each thread copies WORK_PER_THREAD values to the local
        // memory
        for (int w = 0; w < WORK_PER_THREAD; w++) {
          const int tiled_i = THREAD_BLOCK_SIZE * tile_ind + thread_block_row;
          const int tiled_j = THREAD_BLOCK_SIZE * tile_ind + thread_block_col;
          A_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                 [thread_block_row]
              = C[t * temp_rows * temp_rows
                  + (tiled_j + w * THREAD_BLOCK_SIZE_COL) * temp_rows + i];
          if ((offset + j + w * THREAD_BLOCK_SIZE_COL) <= (tiled_i
                    + offset)) {
            B_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                  [thread_block_row]
                = A[(offset + j + w * THREAD_BLOCK_SIZE_COL) * M + tiled_i
                    + offset];
          } else {
            B_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                  [thread_block_row]
                = 0;
          }
        }
        // wait until all tile values are loaded to the local memory
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int block_ind = 0; block_ind < THREAD_BLOCK_SIZE; block_ind++) {
          for (int w = 0; w < WORK_PER_THREAD; w++) {
            acc[w] += A_local[block_ind][thread_block_row]
                      * B_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                               [block_ind];
          }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      // save the values
      for (int w = 0; w < WORK_PER_THREAD; w++) {
        // each thread saves WORK_PER_THREAD values
        if ((i + temp_rows + offset) < non_padded_rows) {
          A[(offset + j + w * THREAD_BLOCK_SIZE_COL) * M
              + i + temp_rows + offset] = -acc[w];
        }
      }
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/matrix_multiply.hpp add() \endlink
 */
const local_range_kernel<cl::Buffer, cl::Buffer, int, int, int, int>
    lower_tri_inverse_step3("lower_tri_inverse_step3",
                            lower_tri_inverse_step3_kernel_code,
                            {{"THREAD_BLOCK_SIZE", 32},
                             {"WORK_PER_THREAD", 8}});
}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
