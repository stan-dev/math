#ifndef STAN_MATH_GPU_KERNELS_EIGENDECOMPOSITION_HPP
#define STAN_MATH_GPU_KERNELS_EIGENDECOMPOSITION_HPP
#ifdef STAN_OPENCL

#include <stan/math/gpu/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
const char* eigendecomp_householder_kernel_code = STRINGIFY(
// \endcond
        /**
         * Calculates householder vector.
         * Must be run with 1 workgroup of 128 threads.
         */
        __kernel void eigendecomp_householder(const __global double* P, __global double* U, __global double* V,
                                              const int P_start, const int P_end, const int j, const int UV_rows) {
          const int lid = get_local_id(0);
          const int gid = get_global_id(0);
          const int gsize = get_global_size(0);
          const int lsize = get_local_size(0);
          const int ngroups = get_num_groups(0);
          const int wgid = get_group_id(0);

          double q = 0;

          const int P_span = P_end - P_start;
          for (int i = lid; i < P_span; i += lsize) {
            double acc = 0;
            if (j != 0) {
              for (int k = 0; k < j; k++) {
                acc += U[UV_rows * k + j - 1 + i] * V[UV_rows * k + j - 1] +
                       V[UV_rows * k + j - 1 + i] * U[UV_rows * k + j - 1];
              }
            }
            double tmp = P[P_start + i] - acc;
            P[P_start + i] = tmp;
            q += tmp * tmp;
          }
          __local double norm_local[128];
          q_local[lid] = q;
          barrier(CLK_LOCAL_MEM_FENCE);
          if (gid == 0) {
            for (int i = 1; i < 128; i++) {
              q += q_local[i];
            }

            double p1 = P[P_start + 1];
            double alpha = -copysign(sqrt(q), p1);
            q -= p1 * p1;
            p1 -= alpha;
            P[P_start + 1] = p1;
            q += p1 * p1;
            q_local[0] = sqrt(q);
            q_local[1] = alpha;
            V[UV_rows*j]=q/sqrt(2.);
          }
          barrier(CLK_LOCAL_MEM_FENCE);
          q = q_local[0];
          alpha = q_local[1];
          double multi = sqrt(2.) / q;
          for (int i = lid + 1; i < P_span; i += lsize) {
            P[i] *= multi;
          }
        }
// \cond
);
// \endcond

// \cond
const char* eigendecomp_apply_Q1_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Calculates res = Pb^T * Pc, where Pb is a block of matrix P, and Pc is a column of matrix P. Both only from start_row down.
     * Launch 1 work group per element of result vector. Local size must be 64.
     * @param P Matrix P.
     * @param[out] res The result
     * @param total_rows Number of rows of matrix P.
     * @param total_cols Number of columns of matrix P.
     * @param start_row First row of the matrix to use.
     * @param start_col First column of the block to use for left side of the multiplication.
     * @param end_col Column to use for right side of multiplication, but not left side.
     */
    __kernel void eigendecomp_apply_Q1(const __global double* P, __global double* res, const int total_rows,
            const int total_cols, const int start_row, const int start_col, const int end_col) {
      const int lid = get_local_id(0);
      const int gid = get_global_id(0);
      const int gsize = get_global_size(0);
      const int lsize = get_local_size(0);
      const int ngroups = get_num_groups(0);
      const int wgid = get_group_id(0);
      __local double res_loc[64];
      double acc = 0;
      int ncols = end_col - start_col;

      const __global double* M = P + total_rows * (start_col + wgid);
      const __global double* V = P + total_rows * end_col;
      for (int i = start_row + lid; i < total_rows; i += 64) { //go over column of the matrix in steps of 64
        acc += M[i] * V[i];
        //printf("%d using %lf, \t %lf\n",gid,M[i], V[i]);
      }
      res_loc[lid]=acc;
      barrier(CLK_LOCAL_MEM_FENCE);
      if(lid==0){
        for(int i=1;i<lsize;i++){
          acc+=res_loc[i];
        }
        res[wgid]=acc;
      }
    }
    // \cond
);
// \endcond

// \cond
const char* eigendecomp_apply_Q1v2_kernel_code = STRINGIFY(
// \endcond
/**
 * Calculates res = Pb^T * Pc, where Pb is a block of matrix P, and Pc is a column of matrix P. Both only from start_row down.
 * Local size must be <=64.
 * @param P Matrix P.
 * @param[out] res The result
 * @param total_rows Number of rows of matrix P.
 * @param total_cols Number of columns of matrix P.
 * @param start_row First row of the matrix to use.
 * @param start_col First column of the block to use for left side of the multiplication.
 * @param end_col Column to use for right side of multiplication, but not left side.
 */
        __kernel void eigendecomp_apply_Q1(const __global double* P, __global double* res, const int total_rows,
                                           const int total_cols, const int start_row, const int start_col, const int end_col) {
          const int lid = get_local_id(0);
          const int gid = get_global_id(0);
          const int gsize = get_global_size(0);
          const int lsize = get_local_size(0);
          const int ngroups = get_num_groups(0);
          const int wgid = get_group_id(0);
          __local double res_loc[64];
          int ncols = end_col - start_col;

          const __global double* V = P + total_rows * end_col;
          const __global double* M = P;
          for(int col=start_col + wgid; col<end_col; col+=ngroups) {
            double acc = 0;
            for (int i = start_row + lid; i < total_rows; i += lsize) { //go over column of the matrix in steps of 64
              acc += M[total_rows *col+i] * V[i];
              //printf("%d using %lf, \t %lf\n",gid,M[i], V[i]);
            }
            res_loc[lid] = acc;
            barrier(CLK_LOCAL_MEM_FENCE);
            if (lid == 0) {
              for (int i = 1; i < lsize; i++) {
                acc += res_loc[i];
              }
              res[col] = acc;
            }
            barrier(CLK_LOCAL_MEM_FENCE);
          }
        }
// \cond
);
// \endcond

// \cond
const char* eigendecomp_apply_Q2_kernel_code = STRINGIFY(
// \endcond
    /**
     * Calculates Mc = Vc - Ml * v, where Mc is the ncols-th column of the matrix M, Vc is the ncols-th column of the matrix V,
     * Ml is left part (ncols columns) of the matrix M and v is a vector. Instead of first ncols elements of Vc zeros are used.
     * Launch 1 thread per element of result vector, rounded up to multiple of 64. Local size must be 64.
     * @param[in,out] M Matrix M.
     * @param v Vector v.
     * @param V Matrix V.
     * @param total_rows Number of rows of M and v and V.
     * @param total_cols Number of columns of M.
     * @param ncols Number of columns of M to use for multiplication.
     */
    __kernel void eigendecomp_apply_Q2(__global double* M, const __global double* v, const __global double* V,
            const int total_rows, const int total_cols, const int ncols, const int k) {
      const int lid = get_local_id(0);
      const int gid = get_global_id(0);
      const int gsize = get_global_size(0);
      const int lsize = get_local_size(0);
      const int ngroups = get_num_groups(0);
      const int wgid = get_group_id(0);
      __local double v_loc[64];
      double acc = 0;

      M += gid;
      __global double* M_start=M;
      v += lid;
      for (int k = 0; k < ncols; k += 64) {
        int end = min(64,ncols-k);
        if(lid<end){
          v_loc[lid] = *v;
          v += 64;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        if(gid < total_rows){
          for (int j = 0; j < end; j++) {
            //printf("%d using %lf, \t %lf\n",gid,*M, v_loc[j]);
            acc += *M * v_loc[j];
            M += total_rows;
          }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      if(gid<ncols){
        //printf("%d acc: %lf\n",gid, acc);
        M_start[total_rows*ncols] = - acc;
      }
      else if(gid<total_rows){
        int V_index = (total_rows + k + 1) * (ncols + k) + k + 1 + gid;
        //printf("%d acc: %lf,  V[]: %lf\n", gid, acc, V[V_index]);
        M_start[total_rows*ncols] = V[V_index] - acc;
      }
    }
// \cond
);
// \endcond

// \cond
const char* subtract_matrix_multiply_kernel_code = STRINGIFY(
// \endcond
/**
 * Calculates Cb -= A * B, where Cb is bottom part of the matrix C, while A and B are matrices.
 *
 * @param[in] A the left matrix in matrix multiplication
 * @param[in] B the right matrix in matrix multiplication
 * @param[out] C the output matrix
 * @param[in] M Number of rows for matrix A and bottom rows for C
 * @param[in] N Number of cols for matrix B
 * @param[in] K Number of cols for matrix A and number of rows for matrix B
 * @param[int] toatal_rows Number of rows for matrix C
 */
        __kernel void matrix_multiply(const __global double* A,
                                      const __global double* B, __global double* C,
                                      const int M, const int N, const int K, const int total_rows) {
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
                          * B_local[thread_block_col + w * THREAD_BLOCK_SIZE_COL]
                          [block_ind];
              }
            }
            barrier(CLK_LOCAL_MEM_FENCE);
          }
          int skip_rows = total_rows - M;
          // save the values
          for (int w = 0; w < WORK_PER_THREAD; w++) {
            // each thread saves WORK_PER_THREAD values
            C[(j + w * THREAD_BLOCK_SIZE_COL) * total_rows + skip_rows + i] = acc[w];
          }
        }
// \cond
);
// \endcond

/**
 * See the docs for \link kernels/matrix_multiply.hpp add() \endlink
 */
//const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, int, int, int>
//    matrix_multiply("matrix_multiply", matrix_multiply_kernel_code,
//                    {{"THREAD_BLOCK_SIZE", 32}, {"WORK_PER_THREAD", 8}});

const local_range_kernel<cl::Buffer, cl::Buffer, int, int, int, int, int>
    eigendecomp_apply_Q1("eigendecomp_apply_Q1", eigendecomp_apply_Q1v2_kernel_code);

const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, int, int, int, int>
        eigendecomp_apply_Q2("eigendecomp_apply_Q2", eigendecomp_apply_Q2_kernel_code);

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
