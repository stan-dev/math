#ifndef STAN_MATH_GPU_KERNELS_TRIDIAGONALIZATION_HPP
#define STAN_MATH_GPU_KERNELS_TRIDIAGONALIZATION_HPP

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
          if(q != 0) {
            double multi = sqrt(2.) / q;
            for (int i = lid + 1; i < P_span; i += lsize) {
              P[P_start + i] *= multi;
            }
          }
          if(gid==0){
            P[P_rows * (k + j + 1) + k + j] = P[P_rows * (k + j) + k + j + 1] * q / sqrt(2.) + alpha;
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

const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, int, int, int>
        tridiagonalization_v1("tridiagonalization_v1", tridiagonalization_v1_kernel_code);

const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, int, int, int, int>
        tridiagonalization_v2("tridiagonalization_v2", tridiagonalization_v2_kernel_code);

const local_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, int, int, int, int>
        tridiagonalization_v3("tridiagonalization_v3", tridiagonalization_v3_kernel_code);

const global_range_kernel<cl::Buffer, cl::Buffer, int, int, int>
        subtract_twice("subtract_twice", subtract_twice_kernel_code);

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
