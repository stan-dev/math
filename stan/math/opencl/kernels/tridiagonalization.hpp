#ifndef STAN_MATH_GPU_KERNELS_TRIDIAGONALIZATION_HPP
#define STAN_MATH_GPU_KERNELS_TRIDIAGONALIZATION_HPP

#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static const char* tridiagonalization_householder_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Calculates householder vector and first element of the vector v.
     * Must be run with 1 workgroup of LOCAL_SIZE_ threads.
     * @param[in,out] P packed matrix being constructed
     * @param[in,out] V matrix V
     * @param[out] q_glob q
     * @param P_rows Number of rows of the packed matrix
     * @param V_rows Number of rows of the matrix V
     * @param j Start column of the block to work on
     * @param k Index of the householder vector in the block to create
     */
    __kernel void tridiagonalization_householder(
        __global double* P, __global double* V, __global double* q_glob,
        const int P_rows, const int V_rows, const int j, const int k) {
      const int lid = get_local_id(0);
      const int gid = get_global_id(0);
      const int gsize = get_global_size(0);
      const int lsize = get_local_size(0);
      const int ngroups = get_num_groups(0);
      const int wgid = get_group_id(0);

      double q = 0;

      const int P_start = P_rows * (k + j) + k + j;
      const int P_span = P_rows * (k + j + 1) - P_start;
      for (int i = lid; i < P_span; i += lsize) {
        double acc = 0;
        // apply previous householder reflections from current block to the
        // column we are making the Householder vector from
        for (int l = 0; l < j; l++) {
          acc += P[P_rows * (k + l) + k + j + i] * V[V_rows * l + j - 1]
                 + V[V_rows * l + j - 1 + i] * P[P_rows * (k + l) + k + j];
        }
        double tmp = P[P_start + i] - acc;
        P[P_start + i] = tmp;
        if (i != 0) {
          q += tmp * tmp;
        }
      }
      // calculate column norm between threads
      __local double q_local[LOCAL_SIZE_];
      q_local[lid] = q;
      barrier(CLK_LOCAL_MEM_FENCE);
      for (int step = lsize / REDUCTION_STEP_SIZE; step > 0;
           step /= REDUCTION_STEP_SIZE) {
        if (lid < step) {
          for (int i = 1; i < REDUCTION_STEP_SIZE; i++) {
            q_local[lid] += q_local[lid + step * i];
          }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }

      double alpha;
      if (lid == 0) {
        q = q_local[0];
        double p1 = P[P_start + 1];
        // make Householder vector
        alpha = -copysign(sqrt(q), P[P_start]);
        q -= p1 * p1;
        p1 -= alpha;
        P[P_start + 1] = p1;
        q += p1 * p1;
        q = sqrt(q);
        q_local[0] = q;
        q_local[1] = alpha;
        *q_glob = q;
      }
      barrier(CLK_LOCAL_MEM_FENCE);
      q = q_local[0];
      alpha = q_local[1];
      if (q != 0) {
        double multi = M_SQRT2 / q;
        // normalize the Householder vector
        for (int i = lid + 1; i < P_span; i += lsize) {
          P[P_start + i] *= multi;
        }
      }
      if (gid == 0) {
        P[P_rows * (k + j + 1) + k + j]
            = P[P_rows * (k + j) + k + j + 1] * q / M_SQRT2 + alpha;
      }
    }
    // \cond
);
// \endcond

// \cond
static const char* tridiagonalization_v_step_1_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Calculates first part of constructing the vector v: Uu = Pb * u and Vu =
     * Vl * u. Pb is a block of packed matrix, Vl is left part of matrix V and u
     * is householder vector. Must be run with number of work groups equal to
     * size of resulting vectors and 64 threads per work group.
     * @param P Packed matrix being constructed.
     * @param V Matrix V.
     * @param[out] Uu First resulting vector.
     * @param[out] Vu Second resulting vector.
     * @param P_rows Number of rows of the packed matrix
     * @param V_rows Number of rows of the matrix V
     * @param k Index of the householder vector in the block we use as input
     */
    __kernel void tridiagonalization_v_step_1(
        const __global double* P, const __global double* V, __global double* Uu,
        __global double* Vu, const int P_rows, const int V_rows, const int k) {
      const int lid = get_local_id(0);
      const int gid = get_global_id(0);
      const int gsize = get_global_size(0);
      const int lsize = get_local_size(0);
      const int ngroups = get_num_groups(0);
      const int wgid = get_group_id(0);

      __local double res_loc1[LOCAL_SIZE_];
      __local double res_loc2[LOCAL_SIZE_];
      double acc1 = 0;
      double acc2 = 0;

      const __global double* vec = P + P_rows * (k + ngroups) + k + ngroups + 1;
      const __global double* M1 = P + P_rows * (k + wgid) + k + ngroups + 1;
      const __global double* M2 = V + V_rows * wgid + ngroups;
      for (int i = lid; i < P_rows - k - ngroups - 1;
           i += LOCAL_SIZE_) {  // go over column of the matrix in steps of 64
        double v = vec[i];
        acc1 += M1[i] * v;
        acc2 += M2[i] * v;
      }
      res_loc1[lid] = acc1;
      res_loc2[lid] = acc2;
      barrier(CLK_LOCAL_MEM_FENCE);

      for (int step = lsize / REDUCTION_STEP_SIZE; step > 0;
           step /= REDUCTION_STEP_SIZE) {
        if (lid < step) {
          for (int i = 1; i < REDUCTION_STEP_SIZE; i++) {
            res_loc1[lid] += res_loc1[lid + step * i];
            res_loc2[lid] += res_loc2[lid + step * i];
          }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      if (lid == 0) {
        Uu[wgid] = res_loc1[0];
        Vu[wgid] = res_loc2[0];
      }
    }
    // \cond
);
// \endcond

// \cond
static const char* tridiagonalization_v_step_2_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Second part in constructing vector v: v = Pb * u + V * Uu + U * Vu. Pb is
     * a block of packed matrix and U is householder vector. Pb is symmetric
     * with only lower triangle having values. That is why two columns of V are
     * written that must be added to obtain the vector v. Must be run with 64
     * threads per work group and total number of threads equal or greater than
     * size of result vector.
     * @param P Packed matrix being constructed.
     * @param V Matrix V.
     * @param Uu Uu from previous kernel.
     * @param Vu Vu from previous kernel.
     * @param P_rows Number of rows of the packed matrix
     * @param V_rows Number of rows of the matrix V
     * @param k Index of the householder vector in the block we use as input
     * @param j Start column of the block to work on
     */
    __kernel void tridiagonalization_v_step_2(
        const __global double* P, __global double* V, const __global double* Uu,
        const __global double* Vu, const int P_rows, const int V_rows,
        const int k, const int j) {
      const int lid = get_local_id(0);
      const int gid = get_global_id(0);
      const int gsize = get_global_size(0);
      const int lsize = get_local_size(0);
      const int ngroups = get_num_groups(0);
      const int wgid = get_group_id(0);

      int work = P_rows - k - j - 1;
      double acc = 0;

      const __global double* vec = P + P_rows * (k + j) + k + j + 1;
      const __global double* M1 = P + P_rows * (k + j + 1) + k + j + 1;
      const __global double* M2 = P + P_rows * k + k + j + 1;
      const __global double* M3 = V + j;
      int i;
      if (gid < work) {
        for (i = 0; i <= gid; i++) {
          acc += M1[P_rows * i + gid] * vec[i];
        }
        for (int i = 0; i < j; i++) {
          acc -= M2[P_rows * i + gid] * Vu[i];
          acc -= M3[V_rows * i + gid] * Uu[i];
        }
        V[V_rows * j + gid + j] = acc;
      }
      float work_per_group
          = (float)work / ngroups;  // NOLINT(readability/casting)
      int start = work_per_group * wgid;
      int end = work_per_group * (wgid + 1);
      __local double res_loc[LOCAL_SIZE_];
      for (int i = start; i < end; i += 1) {
        acc = 0;
        for (int l = i + 1 + lid; l < work; l += LOCAL_SIZE_) {
          acc += M1[P_rows * i + l] * vec[l];
        }
        res_loc[lid] = acc;
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int step = lsize / REDUCTION_STEP_SIZE; step > 0;
             step /= REDUCTION_STEP_SIZE) {
          if (lid < step) {
            for (int i = 1; i < REDUCTION_STEP_SIZE; i++) {
              res_loc[lid] += res_loc[lid + step * i];
            }
          }
          barrier(CLK_LOCAL_MEM_FENCE);
        }
        if (lid == 0) {
          V[V_rows * (j + 1) + i + j] = res_loc[lid];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
    }
    // \cond
);
// \endcond

// \cond
static const char* tridiagonalization_v_step_3_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Third part in constructing vector v: v-=0.5*(v^T*u)*u, where u is the
     * householder vector.
     * @param[in,out] P packed matrix being constructed
     * @param[in,out] V matrix V
     * @param[out] q q
     * @param P_rows Number of rows of the packed matrix
     * @param V_rows Number of rows of the matrix V
     * @param k Index of the householder vector in the block to create
     * @param j Start column of the block to work on
     */
    __kernel void tridiagonalization_v_step_3(
        __global double* P, __global double* V, __global double* q,
        const int P_rows, const int V_rows, const int k, const int j) {
      const int lid = get_local_id(0);
      const int gid = get_global_id(0);
      const int gsize = get_global_size(0);
      const int lsize = get_local_size(0);
      const int ngroups = get_num_groups(0);
      const int wgid = get_group_id(0);

      __global double* u = P + P_rows * (k + j) + k + j + 1;
      __global double* v = V + V_rows * j + j;
      double acc = 0;

      for (int i = lid; i < P_rows - k - j - 1; i += LOCAL_SIZE_) {
        double vi = v[i] + v[i + V_rows];
        v[i] = vi;
        acc += u[i] * vi;
      }
      __local double res_loc[LOCAL_SIZE_];
      res_loc[lid] = acc;
      barrier(CLK_LOCAL_MEM_FENCE);
      for (int step = lsize / REDUCTION_STEP_SIZE; step > 0;
           step /= REDUCTION_STEP_SIZE) {
        if (lid < step) {
          for (int i = 1; i < REDUCTION_STEP_SIZE; i++) {
            res_loc[lid] += res_loc[lid + step * i];
          }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      acc = res_loc[0] * 0.5;
      for (int i = lid; i < P_rows - k - j - 1; i += LOCAL_SIZE_) {
        v[i] -= acc * u[i];
      }
      if (gid == 0) {
        P[P_rows * (k + j + 1) + k + j] -= *q / M_SQRT2 * u[0];
      }
    }
    // \cond
);
// \endcond

const kernel_cl<in_out_buffer, in_out_buffer, out_buffer, int, int, int, int>
    tridiagonalization_householder("tridiagonalization_householder",
                                   {tridiagonalization_householder_kernel_code},
                                   {{"REDUCTION_STEP_SIZE", 4},
                                    {"LOCAL_SIZE_", 1024}});

const kernel_cl<in_buffer, in_buffer, out_buffer, out_buffer, int, int, int>
    tridiagonalization_v_step_1("tridiagonalization_v_step_1",
                                {tridiagonalization_v_step_1_kernel_code},
                                {{"REDUCTION_STEP_SIZE", 4},
                                 {"LOCAL_SIZE_", 64}});

const kernel_cl<in_buffer, out_buffer, in_buffer, in_buffer, int, int, int, int>
    tridiagonalization_v_step_2("tridiagonalization_v_step_2",
                                {tridiagonalization_v_step_2_kernel_code},
                                {{"REDUCTION_STEP_SIZE", 4},
                                 {"LOCAL_SIZE_", 64}});

const kernel_cl<in_out_buffer, in_out_buffer, out_buffer, int, int, int, int>
    tridiagonalization_v_step_3("tridiagonalization_v_step_3",
                                {tridiagonalization_v_step_3_kernel_code},
                                {{"REDUCTION_STEP_SIZE", 4},
                                 {"LOCAL_SIZE_", 1024}});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
