#ifndef STAN_MATH_OPENCL_KERNELS_MATRIX_MULTIPLY_HPP
#define STAN_MATH_OPENCL_KERNELS_MATRIX_MULTIPLY_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static const char* matrix_multiply_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Matrix multiplication on the OpenCL device
     *
     * @param[in] A the left matrix in matrix multiplication
     * @param[in] B the right matrix in matrix multiplication
     * @param[out] C the output matrix
     * @param[in] M Number of rows for matrix A
     * @param[in] N Number of cols for matrix B
     * @param[in] K Number of cols for matrix A and number of rows for matrix B
     * @param[in] lower_upper_A the triangularity of A (lower, upper or none)
     * @param[in] lower_upper_B the triangularity of B (lower, upper or none)
     */
    __kernel void matrix_multiply(
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

      // This kernel is based on the well known
      // general matrix multiplication kernels that
      // use tiling for shared memory
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
      // is then a general matrix multiply
      const int end_tile_A
          = lower_upper_A == LOWER ? (i / THREAD_BLOCK_SIZE) : (num_tiles - 1);
      const int end_tile_B
          = lower_upper_B == UPPER ? (j / THREAD_BLOCK_SIZE) : (num_tiles - 1);
      const int start_tile_A
          = lower_upper_A == UPPER ? (i / THREAD_BLOCK_SIZE) : 0;
      const int start_tile_B
          = lower_upper_B == LOWER ? (j / THREAD_BLOCK_SIZE) : 0;
      const int start_tile = max(start_tile_A, start_tile_B);
      const int end_tile = min(end_tile_A, end_tile_B);  // NOLINT

      for (int tile_idx = start_tile; tile_idx <= end_tile; tile_idx++) {
        const int tiled_i = THREAD_BLOCK_SIZE * tile_idx + thread_block_row;
        const int tiled_j = THREAD_BLOCK_SIZE * tile_idx + thread_block_col;
        // each thread copies WORK_PER_THREAD values to the local
        // memory
        for (int w = 0; w < WORK_PER_THREAD; w++) {
          // For the tiles on the diagonal we can ignore the values over
          // the diagonal if the matrix is lower triangular or under
          // the diagonal if the matrix is upper triangular
          A_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                 [thread_block_row]
              = A[(tiled_j + w * THREAD_BLOCK_SIZE_COL) * M + i];
          B_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                 [thread_block_row]
              = B[(j + w * THREAD_BLOCK_SIZE_COL) * K + tiled_i];
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
 * See the docs for \link kernels/matrix_multiply.hpp matrix_multiply() \endlink
 */
const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, int, int, int,
                         TriangularViewCL, TriangularViewCL>
    matrix_multiply("matrix_multiply", matrix_multiply_kernel_code,
                    {{"THREAD_BLOCK_SIZE", 32}, {"WORK_PER_THREAD", 8}});

// \cond
static const char* matrix_vector_multiply_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Matrix-vector multiplication R=A*B on the OpenCL device
     *
     * @param[in] A matrix in matrix-vector multiplication
     * @param[in] B vector in matrix-vector multiplication
     * @param[out] R the output vector
     * @param[in] M Number of rows for matrix A
     * @param[in] N Number of cols for matrix A and number of rows for vector B
     */
    __kernel void matrix_vector_multiply(
        const __global double* A, const __global double* B, __global double* R,
        const int M, const int N) {
      const int gid = get_global_id(0);

      double acc = 0;
      for (int i = 0, j = 0; i < N; i++, j += M) {
        acc += A[j + gid] * B[i];
      }
      R[gid] = acc;
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/matrix_multiply.hpp matrix_vector_multiply()
 * \endlink
 */
const global_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, int, int>
    matrix_vector_multiply("matrix_vector_multiply",
                           matrix_vector_multiply_kernel_code);

// \cond
static const char* row_vector_matrix_multiply_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Row vector-matrix multiplication R=A*B on the OpenCL device
     *
     * @param[in] A row vector in row vector-matrix multiplication
     * @param[in] B matrix in row vector-matrix multiplication
     * @param[out] R the output vector
     * @param[in] N Number of cols for vector A and number of rows for matrix B
     * @param[in] K Number of cols for matrix B
     */
    __kernel void row_vector_matrix_multiply(
        const __global double* A, const __global double* B, __global double* R,
        const int N, const int K) {
      const int lid = get_local_id(0);
      const int gid = get_global_id(0);
      const int wgid = get_group_id(0);

      double acc = 0;
      for (int i = lid; i < N; i += 64) {
        acc += A[i] * B[i + wgid * N];
      }
      __local double res_loc[64];
      res_loc[lid] = acc;
      barrier(CLK_LOCAL_MEM_FENCE);
      if (lid < 16) {
        res_loc[lid]
            += res_loc[lid + 16] + res_loc[lid + 32] + res_loc[lid + 48];
      }
      barrier(CLK_LOCAL_MEM_FENCE);
      if (lid < 4) {
        res_loc[lid] += res_loc[lid + 4] + res_loc[lid + 8] + res_loc[lid + 12];
      }
      barrier(CLK_LOCAL_MEM_FENCE);
      if (lid == 0) {
        acc = res_loc[lid] + res_loc[lid + 1] + res_loc[lid + 2]
              + res_loc[lid + 3];
        R[wgid] = acc;
      }
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/matrix_multiply.hpp
 * row_vector_matrix_multiply() \endlink
 */
const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, int, int>
    row_vector_matrix_multiply("row_vector_matrix_multiply",
                               row_vector_matrix_multiply_kernel_code);

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
