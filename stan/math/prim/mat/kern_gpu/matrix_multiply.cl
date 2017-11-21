R"=====(
#ifdef cl_khr_fp64
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
    #error "Double not supported by OpenCL implementation."
#endif

#define a(i,j)  a[j*rows+i]
#define b(i,j)  b[j*rows+i]
#define c(i,j)  c[j*rows+i]

__kernel void scalar_mul_diagonal(
          __global double *a,
          double scalar,
          unsigned int rows,
          unsigned int cols) {
    int i = get_global_id(0);
    if ( i < rows && i < cols ) {
     a(i,i) *= scalar;
    }
}

__kernel void scalar_mul(
          __global double *a,
          __global double *b,
          double scalar,
          unsigned int rows,
          unsigned int cols) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    if ( i < rows && j < cols ) {
     a(i,j) = b(i,j)*scalar;
    }
}

#define A(i,j)  A[j*M+i]
#define B(i,j)  B[j*K+i]
#define C(i,j)  C[j*M+i]

#define TS 16
__kernel void basic_multiply(const int M, const int N, const int K,
                      const __global double* A,
                      const __global double* B,
                      __global double* C) {
        
    const int local_i = get_local_id(0);
    const int local_j = get_local_id(1);
    const int i = TS*get_group_id(0) + local_i;
    const int j = TS*get_group_id(1) + local_j;

    __local double Asub[TS][TS];
    __local double Bsub[TS][TS];

    double acc = 0.0;

    const int numTiles = (K+TS-1)/TS;
    
    for (int t=0; t < numTiles; t++) {
        
        const int tiled_i = TS*t + local_i;
        const int tiled_j = TS*t + local_j;

        if ( i < M && tiled_j < K ) {
         Asub[local_j][local_i] = A(i, tiled_j);
        } else {
         Asub[local_j][local_i] = 0.0;
        }

        if ( tiled_i < K && j < N ) {
         Bsub[local_j][local_i] = B(tiled_i, j);
        } else {
         Bsub[local_j][local_i] = 0.0;
        }

        barrier(CLK_LOCAL_MEM_FENCE);

        for (int k=0; k < TS; k++) {
            acc += Asub[k][local_i] * Bsub[local_j][k];
        }

        barrier(CLK_LOCAL_MEM_FENCE);
    }

    if ( j < N && i < M ) {
     C(i,j) = acc;
    }
}
 
#define TS11 16
__kernel void multiply_self_transposed(const int M, const int N, const int K,
                      const __global double* A,
                      const __global double* B,
                      __global double* C) {
    
    const int local_i = get_local_id(0);
    const int local_j = get_local_id(1);
    const int i = TS11*get_group_id(0) + local_i;
    const int j = TS11*get_group_id(1) + local_j;
    const int jMin = TS11*get_group_id(1);
    const int iMax = TS11*get_group_id(0)  + get_local_size(0);

    __local double Asub[TS11][TS11];
    __local double Bsub[TS11][TS11];

    double acc = 0.0;

    const int numTiles = (K+TS11-1)/TS11;
    for (int t=0; t < numTiles; t++) {
        
        const int tiled_i = TS11*t + local_i;
        const int tiled_j = TS11*t + local_j;
        if(jMin <= iMax){
          
          if ( i < M && tiled_j < K ) {
           Asub[local_j][local_i] = A(i, tiled_j);
          } else {
           Asub[local_j][local_i] = 0.0;
          }
          if ( tiled_i < K && j < N ) {
           Bsub[local_j][local_i] = B(tiled_i, j);
          } else {
           Bsub[local_j][local_i] = 0.0;
          }
        }
        
        barrier(CLK_LOCAL_MEM_FENCE);

        
        if( j <= i ){
          for (int k=0; k < TS11; k++) {
              acc += Asub[k][local_i] * Bsub[local_j][k];
          }
        }

        barrier(CLK_LOCAL_MEM_FENCE);
    }

    if ( j < N && i < M && j <= i ) {
     C(i,j) = acc;
     C(j,i) = acc;
    }
 }
)====="
