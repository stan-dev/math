#ifndef STAN_MATH_GPU_KERNELS_LOWER_TRI_INVERSE_STEP3_HPP
#define STAN_MATH_GPU_KERNELS_LOWER_TRI_INVERSE_STEP3_HPP
#ifdef STAN_OPENCL

#include <stan/math/gpu/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
const char* negative_rectangular_lower_triangular_multiply_kernel_code
    = STRINGIFY(
        // \endcond
        /**
         * Calculates C = -B * A where B is rectangular and A is a lower
         * triangular.
         * The full inverse requires calculation of the lower left rectangular
         * matrix within the lower left triangular C3 = -C2*A3*C1. where C2 is the
         * inverse of the bottom right lower triangular, C1 is the inverse of the
         * upper left lower and A3 is the original lower triangulars lower left
         * rectangular. This kernel performs multiplications on submatrices in
         * the input matrix A in parallel and includes optimizations
         * to account for the lower triangular input matrix A.
         *
         *
         * @param[in, out] A Input matrix that is being inverted.
         * @param[in] temp Temporary matrix with the intermediate results.
         * @param A_rows Number of rows for A.
         * @param rows The number of rows in a single matrix of the batch
         * @note Code is a <code>const char*</code> held in
         *  negative_rectangular_lower_triangular_multiply_kernel_code
         *  Used in math/gpu/lower_tri_inverse.hpp.
         *  This kernel uses the helper macros available in helpers.cl.
         */
        __kernel void negative_rectangular_lower_triangular_multiply(
            __global double* A, const __global double* temp, const int A_rows,
            const int rows) {
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

          const int num_tiles
              = (rows + THREAD_BLOCK_SIZE - 1) / THREAD_BLOCK_SIZE;
          // iterate over all tiles
          for (int tile_ind = 0; tile_ind < num_tiles; tile_ind++) {
            // each thread copies WORK_PER_THREAD values to the local
            // memory
            for (int w = 0; w < WORK_PER_THREAD; w++) {
              const int tiled_i
                  = THREAD_BLOCK_SIZE * tile_ind + thread_block_row;
              const int tiled_j
                  = THREAD_BLOCK_SIZE * tile_ind + thread_block_col;
              const int temp_global_col = tiled_j + w * THREAD_BLOCK_SIZE_COL;
              const int C1_global_col = offset + j + w * THREAD_BLOCK_SIZE_COL;
              const int C1_global_row = tiled_i + offset;
              const int local_col
                  = thread_block_col + w * THREAD_BLOCK_SIZE_COL;
              const int local_row = thread_block_row;
              if ((temp_global_col) < rows && i < rows) {
                temp_local[local_col][local_row]
                    = temp[result_matrix_id * rows * rows
                           + temp_global_col * rows + i];
              } else {
                temp_local[local_col][local_row] = 0.0;
              }
              if (C1_global_col <= C1_global_row) {
                C1_local[local_col][local_row]
                    = A[C1_global_col * A_rows + C1_global_row];
              } else {
                C1_local[local_col][local_row] = 0;
              }
            }
            // wait until all tile values are loaded to the local memory
            barrier(CLK_LOCAL_MEM_FENCE);
            for (int block_ind = 0; block_ind < THREAD_BLOCK_SIZE;
                 block_ind++) {
              for (int w = 0; w < WORK_PER_THREAD; w++) {
                const int local_col
                    = thread_block_col + w * THREAD_BLOCK_SIZE_COL;
                const int local_row = thread_block_row;
                acc[w] += temp_local[block_ind][local_row]
                          * C1_local[local_col][block_ind];
              }
            }
            barrier(CLK_LOCAL_MEM_FENCE);
          }
          // save the values
          const int A_global_row = i + rows + offset;
          const int A_global_col_offset = offset + j;
          for (int w = 0; w < WORK_PER_THREAD; w++) {
            const int A_global_col
                = A_global_col_offset + w * THREAD_BLOCK_SIZE_COL;
            // each thread saves WORK_PER_THREAD values
            A[A_global_col * A_rows + i + rows + offset] = -acc[w];
          }
        }
        // \cond
    );
// \endcond

/**
 * See the docs for \link kernels/matrix_multiply.hpp add() \endlink
 */
const local_range_kernel<cl::Buffer, cl::Buffer, int, int>
    negative_rectangular_lower_triangular_multiply(
        "negative_rectangular_lower_triangular_multiply",
        negative_rectangular_lower_triangular_multiply_kernel_code,
        {{"THREAD_BLOCK_SIZE", 32}, {"WORK_PER_THREAD", 8}});
}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
