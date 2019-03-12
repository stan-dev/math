#ifndef STAN_MATH_GPU_KERNELS_TRIDIAGONALIZATION_HPP
#define STAN_MATH_GPU_KERNELS_TRIDIAGONALIZATION_HPP
#ifndef STAN_OPENCL
#error "NO STAN_OPENCL"
#endif
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
const char* tridiagonalization_householder_kernel_code = STRINGIFY(
// \endcond
        /**
         * Calculates householder vector and first element of the vector v.
         * Must be run with 1 workgroup of 1024 threads.
         * @param[in,out] P packed matrix being constructed
         * @param[in,out] V matrix V
         * @param[out] q_glob q
         * @param P_rows Number of rows of the packed matrix
         * @param V_rows Number of rows of the matrix V
         * @param j Start column of the block to work on
         * @param k Index of the householder vector in the block to create
         */
        __kernel void tridiagonalization_householder(__global double* P, __global double* V, __global double* q_glob,
                                              const int P_rows, const int V_rows,
                                              const int j, const int k) {
          const int lid = get_local_id(0);
          const int gid = get_global_id(0);
          const int gsize = get_global_size(0);
          const int lsize = get_local_size(0);
          const int ngroups = get_num_groups(0);
          const int wgid = get_group_id(0);

          double q = 0;

          const int P_start = P_rows * (k+j) + k + j;
          const int P_span = P_rows * (k+j+1) - P_start;
          for (int i = lid; i < P_span; i += lsize) {
            double acc = 0;
            if (j != 0) {
              for (int l = 0; l < j; l++) {
                acc += P[P_rows * (k + l) + k + j + i] * V[V_rows * l + j - 1] +
                       V[V_rows * l + j - 1 + i] * P[P_rows * (k + l) + k + j];
              }
            }
            double tmp = P[P_start + i] - acc;
            P[P_start + i] = tmp;
            if(i!=0) {
              q += tmp * tmp;
            }
          }
          __local double q_local[1024];
          q_local[lid] = q;
          barrier(CLK_LOCAL_MEM_FENCE);
          double alpha;
          //efficient parallel reduction
          if(lid<256) {
            q_local[lid] += q_local[lid + 256] + q_local[lid + 512] + q_local[lid + 768];
          }
          barrier(CLK_LOCAL_MEM_FENCE);
          if(lid<64) {
            q_local[lid] += q_local[lid + 64] + q_local[lid + 128] + q_local[lid + 192];
          }
          barrier(CLK_LOCAL_MEM_FENCE);
          if(lid<16) {
            q_local[lid] += q_local[lid + 16] + q_local[lid + 32] + q_local[lid + 48];
          }
          barrier(CLK_LOCAL_MEM_FENCE);
          if(lid<4) {
            q_local[lid] += q_local[lid + 4] + q_local[lid + 8] + q_local[lid + 12];
          }
          barrier(CLK_LOCAL_MEM_FENCE);
          if(lid==0){
            q=q_local[lid]+q_local[lid+1]+q_local[lid+2]+q_local[lid+3];
            double p1 = P[P_start + 1];
            alpha = -copysign(sqrt(q), P[P_start]);
            q -= p1 * p1;
            p1 -= alpha;
            P[P_start + 1] = p1;
            q += p1 * p1;
            q=sqrt(q);
            q_local[0] = q;
            q_local[1] = alpha;
            *q_glob=q;
          }
          barrier(CLK_LOCAL_MEM_FENCE);
          q = q_local[0];
          alpha = q_local[1];
          double multi = sqrt(2.) / q;
          for (int i = lid + 1; i < P_span; i += lsize) {
            P[P_start + i] *= multi;
          }
          if(gid==0){
            P[P_rows * (k + j + 1) + k + j] = P[P_rows * (k + j) + k + j + 1] * q / sqrt(2.) + alpha;
          }
        }
// \cond
);
// \endcond

// \cond
const char* tridiagonalization_householder_v1_kernel_code = STRINGIFY(
// \endcond
        /**
         * Calculates householder vector and first element of the vector v.
         * Must be run with 1 workgroup of 128 threads.
         */
        __kernel void tridiagonalization_householder_v1(__global double* P, __global double* V, __global double* q_glob, __global double* R1, __global double* R2,
                                              const int P_rows, const int V_rows,
                                              const int j, const int k) {
          const int lid = get_local_id(0);
          const int gid = get_global_id(0);
          const int gsize = get_global_size(0);
          const int lsize = get_local_size(0);
          const int ngroups = get_num_groups(0);
          const int wgid = get_group_id(0);

          double q = 0;

          const int P_start = P_rows * (k+j) + k + j;
          const int P_span = P_rows * (k+j+1) - P_start;
          for (int i = lid; i < P_span; i += lsize) {
            double acc = 0;
            if (j != 0) {
              for (int l = 0; l < j; l++) {
                acc += P[P_rows * (k + l) + k + j + i] * V[V_rows * l + j - 1] +
                       V[V_rows * l + j - 1 + i] * P[P_rows * (k + l) + k + j];
              }
            }
            double tmp = P[P_start + i] - acc;
            P[P_start + i] = tmp;
            if(i!=0) {
              q += tmp * tmp;
            }
          }
          __local double q_local[128];
          q_local[lid] = q;
          barrier(CLK_LOCAL_MEM_FENCE);
          double alpha;
          if (gid == 0) {
            for (int i = 1; i < 128; i++) {
              q += q_local[i];
            }

            double p1 = P[P_start + 1];
            alpha = -copysign(sqrt(q), P[P_start]);
            q -= p1 * p1;
            p1 -= alpha;
            P[P_start + 1] = p1;
            q += p1 * p1;
            q=sqrt(q);
            q_local[0] = q;
            q_local[1] = alpha;
            *q_glob=q;
          }
          barrier(CLK_LOCAL_MEM_FENCE);
          q = q_local[0];
          alpha = q_local[1];
          double multi = sqrt(2.) / q;
          for (int i = lid + 1; i < P_span; i += lsize) {
            P[P_start + i] *= multi;
          }
          if(gid==0){
            P[P_rows * (k + j + 1) + k + j] = P[P_rows * (k + j) + k + j + 1] * q / sqrt(2.) + alpha;
          }
          //from here on doing v1 stuff
          barrier(CLK_LOCAL_MEM_FENCE);
          if(j==0){
            return;
          }
          //divide threads in j subgroups - one per result element
          int subgroup_size = 128 / j;
          int working_threads = subgroup_size * j;
          int idx_in_subgroup = gid % subgroup_size;
          int subgroup_idx = gid / subgroup_size;
          const __global double* vec = P + P_rows * (k + j) + k + j + 1;
          const __global double* M1 =  P + P_rows * (k + subgroup_idx) + k + j + 1;
          const __global double* M2 = V + V_rows * subgroup_idx + j;
          double acc1=0;
          double acc2=0;
          if(gid<working_threads) {
            for (int i = idx_in_subgroup; i < P_span - 1; i += subgroup_size) {
              double v = vec[i];
              acc1 += M1[i] * v;
              acc2 += M2[i] * v;
            }
          }
          __local double res_loc1[128];
          __local double res_loc2[128];
          res_loc1[lid]=acc1;
          res_loc2[lid]=acc2;
          barrier(CLK_LOCAL_MEM_FENCE);
          if(gid<working_threads && idx_in_subgroup==0){
            for(int i=lid+1;i<min(lid+subgroup_size, working_threads);i+=1){
              acc1+=res_loc1[i];
              acc2+=res_loc2[i];
            }
            R1[subgroup_idx]=acc1;
            R2[subgroup_idx]=acc2;
          }
        }
// \cond
);
// \endcond

// \cond
const char* tridiagonalization_v1_kernel_code = STRINGIFY(
// \endcond
        /**
         * Calculates first part of constructing the vector v: Uu = Pb * u and Vu = Vl * u. Pb is a block of packed matrix,
         * Vl is left part of matrix V and u is householder vector.
         * Must be run with number of work groups equal to size of resulting vectors and 64 threads per work group.
         * @param P Packed matrix being constructed.
         * @param V Matrix V.
         * @param[out] Uu First resulting vector.
         * @param[out] Vu Second resulting vector.
         * @param P_rows Number of rows of the packed matrix
         * @param V_rows Number of rows of the matrix V
         * @param k Index of the householder vector in the block we use as input
         */
        __kernel void tridiagonalization_v1(const __global double* P, const __global double* V, __global double* Uu, __global double* Vu,
                const int P_rows, const int V_rows, const int k) {
          const int lid = get_local_id(0);
          const int gid = get_global_id(0);
          const int gsize = get_global_size(0);
          const int lsize = get_local_size(0);
          const int ngroups = get_num_groups(0);
          const int wgid = get_group_id(0);

          __local double res_loc1[64];
          __local double res_loc2[64];
          double acc1 = 0;
          double acc2 = 0;

          const __global double* vec = P + P_rows * (k + ngroups) + k + ngroups + 1;
          const __global double* M1 = P + P_rows * (k + wgid) + k + ngroups + 1;
          const __global double* M2 = V + V_rows * wgid + ngroups;
          for (int i = lid; i < P_rows - k - ngroups - 1; i += 64) { //go over column of the matrix in steps of 64
            double v = vec[i];
            acc1 += M1[i] * v;
            acc2 += M2[i] * v;
          }
          res_loc1[lid]=acc1;
          res_loc2[lid]=acc2;
          barrier(CLK_LOCAL_MEM_FENCE);
          if(lid==0){
            for(int i=1;i<lsize;i++){
              acc1+=res_loc1[i];
              acc2+=res_loc2[i];
            }
            Uu[wgid]=acc1;
            Vu[wgid]=acc2;
          }
        }
// \cond
);
// \endcond

// \cond
const char* tridiagonalization_v2_kernel_code = STRINGIFY(
// \endcond
        /**
         * Second part in constructing vector v.
         * @param P Packed matrix being constructed.
         * @param V Matrix V.
         * @param Uu Uu from previous kernel.
         * @param Vu Vu from previous kernel.
         * @param P_rows Number of rows of the packed matrix
         * @param V_rows Number of rows of the matrix V
         * @param k Index of the householder vector in the block we use as input
         * @param j Start column of the block to work on
         */
        __kernel void tridiagonalization_v2(const __global double* P, __global double* V, const __global double* Uu, const __global double* Vu,
                                     const int P_rows, const int V_rows, const int k, const int j) {
          const int lid = get_local_id(0);
          const int gid = get_global_id(0);
          const int gsize = get_global_size(0);
          const int lsize = get_local_size(0);
          const int ngroups = get_num_groups(0);
          const int wgid = get_group_id(0);

          __local double res_loc[64];
          double acc = 0;

          const __global double* vec = P + P_rows * (k + j) + k + j + 1;
          const __global double* M1 = P + P_rows * (k + j + 1) + k + j + 1;
          const __global double* M2 = P + P_rows * k + k + j + 1;
          const __global double* M3 = V + j;
          int i;
          for (i = lid; i < wgid; i += 64) {//non-coalesced
            acc += M1[P_rows * i + wgid] * vec[i];
          }
          for (; i < ngroups; i += 64) {
            acc += M1[P_rows * wgid + i] * vec[i];
          }
          for(int i=lid;i<j;i+=64){
            acc -= M2[P_rows * i + wgid] * Vu[i];
            acc -= M3[V_rows * i + wgid] * Uu[i];
          }
          res_loc[lid]=acc;
          barrier(CLK_LOCAL_MEM_FENCE);
          if(lid==0){
            for(int i=1;i<64;i++){
              acc+=res_loc[i];
            }
            V[V_rows * j + wgid + j]=acc;
          }
        }
// \cond
);
// \endcond

// \cond
const char* tridiagonalization_v2_2_kernel_code = STRINGIFY(
// \endcond
        /**
         * Second part in constructing vector v.
         * @param P Packed matrix being constructed.
         * @param V Matrix V.
         * @param Uu Uu from previous kernel.
         * @param Vu Vu from previous kernel.
         * @param P_rows Number of rows of the packed matrix
         * @param V_rows Number of rows of the matrix V
         * @param k Index of the householder vector in the block we use as input
         * @param j Start column of the block to work on
         */
        __kernel void tridiagonalization_v2(const __global double* P, __global double* V, const __global double* Uu, const __global double* Vu,
                                     const int P_rows, const int V_rows, const int k, const int j) {
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
          if(gid<work) {
            for (i = 0; i <= gid; i++) {
              acc += M1[P_rows * i + gid] * vec[i];
            }
          }
          __local double loc[33][33]; //+1 to avoid bank conflicts
          //handle diagonal blocks
          i=wgid*32;
          if(i+32<work) { //general diagonal case
            for (int l = 0; l < 32; l++) {
              loc[l][lid] = M1[P_rows * (wgid*32 + l) + i + lid];
            }
          }
          else{ //last diagonal block
            for (int l = 0; l < 32; l++) {
              if(wgid*32+l<work && i+lid<work) {
                loc[l][lid] = M1[P_rows * (wgid*32 + l) + i + lid];
              }
            }
          }
          barrier(CLK_LOCAL_MEM_FENCE);

          for(int l=lid+1;l<min(32, work - i);l++){
            acc+=loc[lid][l] * vec[i+l];
          }
          barrier(CLK_LOCAL_MEM_FENCE);

          //general case (blocks going down from diagonal)
          for(i=(wgid+1)*32;i<work-32;i+=32){
            for(int l=0;l<32;l++){
              loc[l][lid] = M1[P_rows * (wgid * 32 + l) + i + lid];
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            for(int l=0;l<32;l++){
              acc+=loc[lid][l] * vec[i + l];
            }
            barrier(CLK_LOCAL_MEM_FENCE);
          }

          //handle last block in column
          if(i<work) { //if not diagonal at same time
            for (int l = 0; l < 32; l++) {
              if (i + lid < work) {
                loc[l][lid] = M1[P_rows * (wgid * 32 + l) + i + lid];
              }
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            for (int l = 0; l < work - i; l++) {
              acc += loc[lid][l] * vec[i + l];
            }
            barrier(CLK_LOCAL_MEM_FENCE);
          }

          if(gid<work) {
            for (int i = 0; i < j; i++) {
              acc -= M2[P_rows * i + gid] * Vu[i];
              acc -= M3[V_rows * i + gid] * Uu[i];
            }
            V[V_rows * j + gid + j] = acc;
          }
        }
// \cond
);
// \endcond

// \cond
const char* tridiagonalization_v2_3_kernel_code = STRINGIFY(
// \endcond
        /**
         * Second part in constructing vector v: v = Pb * u + V * Uu + U * Vu. Pb is a block of packed matrix
         * and U is householder vector. Pb is symmetric with only lower triangel having values. That is why
         * two columns of V are written that must be added to obtain the vector v.
         * Must be run with 64 threads per work group and total number of threads equal or greater than size of
         * result vector.
         * @param P Packed matrix being constructed.
         * @param V Matrix V.
         * @param Uu Uu from previous kernel.
         * @param Vu Vu from previous kernel.
         * @param P_rows Number of rows of the packed matrix
         * @param V_rows Number of rows of the matrix V
         * @param k Index of the householder vector in the block we use as input
         * @param j Start column of the block to work on
         */
        __kernel void tridiagonalization_v2(const __global double* P, __global double* V, const __global double* Uu, const __global double* Vu,
                                            const int P_rows, const int V_rows, const int k, const int j) {
          const int lid = get_local_id(0);
          const int gid = get_global_id(0);
          const int gsize = get_global_size(0);
          const int lsize = get_local_size(0);
          const int ngroups = get_num_groups(0);
          const int wgid = get_group_id(0);

          int work = P_rows - k - j - 1;
          //printf("%d work %d",gid, work);
          double acc = 0;

          const __global double* vec = P + P_rows * (k + j) + k + j + 1;
          const __global double* M1 = P + P_rows * (k + j + 1) + k + j + 1;
          const __global double* M2 = P + P_rows * k + k + j + 1;
          const __global double* M3 = V + j;
          int i;
          if(gid<work) {
            for (i = 0; i <= gid; i++) {
              acc += M1[P_rows * i + gid] * vec[i];
            }
            for (int i = 0; i < j; i++) {
              acc -= M2[P_rows * i + gid] * Vu[i];
              acc -= M3[V_rows * i + gid] * Uu[i];
            }
            V[V_rows * j + gid + j] = acc;
          }
          float work_per_group=(float)work/ngroups;
          int start=work_per_group*wgid;
          int end=work_per_group*(wgid+1);
          for(int i=start;i<end;i+=1){
            __local double res_loc[64];
            acc = 0;
            for(int l=i+1+lid;l<work;l+=64){
              acc += M1[P_rows * i + l] * vec[l];
            }
            res_loc[lid]=acc;
            barrier(CLK_LOCAL_MEM_FENCE);
            if(lid<16) {//efficient parallel reduction
              res_loc[lid] += res_loc[lid + 16] + res_loc[lid + 32] + res_loc[lid + 48];
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            if(lid<4) {
              res_loc[lid] += res_loc[lid + 4] + res_loc[lid + 8] + res_loc[lid + 12];
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            if(lid==0){
              V[V_rows * (j+1) + i + j]=res_loc[lid]+res_loc[lid+1]+res_loc[lid+2]+res_loc[lid+3];
            }
            barrier(CLK_LOCAL_MEM_FENCE);
          }
        }
// \cond
);
// \endcond

// \cond
const char* tridiagonalization_v3_2_kernel_code = STRINGIFY(
// \endcond
        /**
         * Third part in constructing vector v: v-=0.5*(v^T*u)*u, where u is householder vector.
         * @param[in,out] P packed matrix being constructed
         * @param[in,out] V matrix V
         * @param[out] q q
         * @param P_rows Number of rows of the packed matrix
         * @param V_rows Number of rows of the matrix V
         * @param k Index of the householder vector in the block to create
         * @param j Start column of the block to work on
         */
        __kernel void tridiagonalization_v3(__global double* P, __global double* V, __global double* q,
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

          for (int i = lid; i < P_rows - k - j - 1; i += 128) {
            double vi = v[i] + v[i+V_rows];
            v[i]=vi;
            acc+=u[i]*vi;
          }
          __local double res_loc[128];
          res_loc[lid]=acc;
          barrier(CLK_LOCAL_MEM_FENCE);
          if(lid==0){
            for(int i=1;i<128;i++){
              acc+=res_loc[i];
            }
            res_loc[0]=acc;
          }
          barrier(CLK_LOCAL_MEM_FENCE);
          acc=res_loc[0]*0.5;
          for (int i = lid; i < P_rows - k - j - 1; i += 128) {
            v[i] -= acc * u[i];
          }
          if(gid==0) {
            P[P_rows * (k + j + 1) + k + j] -= *q / sqrt(2.) * u[0];
          }
        }
// \cond
);
// \endcond

// \cond
const char* tridiagonalization_v3_kernel_code = STRINGIFY(
// \endcond
        /**
         * Third part in constructing vector v: v-=0.5*(v^T*u)*u, where u is householder vector.
         * @param[in,out] P packed matrix being constructed
         * @param[in,out] V matrix V
         * @param[out] q q
         * @param P_rows Number of rows of the packed matrix
         * @param V_rows Number of rows of the matrix V
         * @param k Index of the householder vector in the block to create
         * @param j Start column of the block to work on
         */
        __kernel void tridiagonalization_v3(__global double* P, __global double* V, __global double* q,
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

          for (int i = lid; i < P_rows - k - j - 1; i += 128) {
            acc+=u[i]*v[i];
          }
          __local double res_loc[128];
          res_loc[lid]=acc;
          barrier(CLK_LOCAL_MEM_FENCE);
          if(lid==0){
            for(int i=1;i<128;i++){
              acc+=res_loc[i];
            }
            res_loc[0]=acc;
          }
          barrier(CLK_LOCAL_MEM_FENCE);
          acc=res_loc[0]*0.5;
          for (int i = lid; i < P_rows - k - j - 1; i += 128) {
            v[i] -= acc * u[i];
          }
          if(gid==0) {
            P[P_rows * (k + j + 1) + k + j] -= *q / sqrt(2.) * u[0];
          }
        }
// \cond
);
// \endcond

// \cond
const char* tridiagonalization_apply_Q1_kernel_code = STRINGIFY(
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
    __kernel void tridiagonalization_apply_Q1(const __global double* P, __global double* res, const int total_rows,
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
const char* tridiagonalization_apply_Q1v2_kernel_code = STRINGIFY(
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
        __kernel void tridiagonalization_apply_Q1(const __global double* P, __global double* res, const int total_rows,
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
const char* tridiagonalization_apply_Q2_kernel_code = STRINGIFY(
// \endcond
    /**
     * Calculates Mc = Vc - Ml * v, where Mc is the ncols-th column of the matrix M, Vc is the ncols-th column of the matrix V,
     * Ml is left part (ncols columns) of the matrix M and v is a vector. Instead of first ncols elements of Vc zeros are used.
     * Must be run with 1 thread per element of result vector, rounded up to multiple of 64. Local size must be 64.
     * @param[in,out] M Matrix M.
     * @param v Vector v.
     * @param V Matrix V.
     * @param total_rows Number of rows of M and v and V.
     * @param total_cols Number of columns of M.
     * @param ncols Number of columns of M to use for multiplication.
     */
    __kernel void tridiagonalization_apply_Q2(__global double* M, const __global double* v, const __global double* V,
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
            acc += *M * v_loc[j];
            M += total_rows;
          }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      if(gid<ncols){
        M_start[total_rows*ncols] = - acc;
      }
      else if(gid<total_rows){
        int V_index = (total_rows + k + 1) * (ncols + k) + k + 1 + gid;
        M_start[total_rows*ncols] = V[V_index] - acc;
      }
    }
// \cond
);
// \endcond

// \cond
const char* subtract_twice_kernel_code = STRINGIFY(
// \endcond
        /**
         * Calculates A -= B + B ^ T, for lower triangular part of bottom right corner of A.
         * @param A[in,out] First matrix.
         * @param B Second matrix.
         * @param A_rows Number of rows of A.
         * @param B_rows Number of rows of B.
         * @param start At which row and column of A to start.
         */
        __kernel void subtract_twice(__global double* A, const __global double* B,
                const int A_rows, const int B_rows, const int start) {
          const int y = get_global_id(0);
          const int x = get_global_id(1);
          if(y>=x){
            A[A_rows * (start + x) + start + y] -= B[B_rows * x + y] + B[B_rows * y + x];
          }
        }
// \cond
);
// \endcond

const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, int, int, int, int>
        tridiagonalization_householder("tridiagonalization_householder", tridiagonalization_householder_kernel_code);

const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, int, int, int, int>
        tridiagonalization_householder_v1("tridiagonalization_householder_v1", tridiagonalization_householder_v1_kernel_code);

const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, int, int, int>
        tridiagonalization_v1("tridiagonalization_v1", tridiagonalization_v1_kernel_code);

const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, int, int, int, int>
        tridiagonalization_v2("tridiagonalization_v2", tridiagonalization_v2_kernel_code);

const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, int, int, int, int>
        tridiagonalization_v2_2("tridiagonalization_v2", tridiagonalization_v2_2_kernel_code);

const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, int, int, int, int>
        tridiagonalization_v2_3("tridiagonalization_v2", tridiagonalization_v2_3_kernel_code);

const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, int, int, int, int>
        tridiagonalization_v3("tridiagonalization_v3", tridiagonalization_v3_kernel_code);

const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, int, int, int, int>
        tridiagonalization_v3_2("tridiagonalization_v3", tridiagonalization_v3_2_kernel_code);

const global_range_kernel<cl::Buffer, cl::Buffer, int, int, int>
        subtract_twice("subtract_twice", subtract_twice_kernel_code);

const local_range_kernel<cl::Buffer, cl::Buffer, int, int, int, int, int>
    tridiagonalization_apply_Q1("tridiagonalization_apply_Q1", tridiagonalization_apply_Q1_kernel_code);

const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, int, int, int, int>
        tridiagonalization_apply_Q2("tridiagonalization_apply_Q2", tridiagonalization_apply_Q2_kernel_code);

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
