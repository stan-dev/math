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

#define WPT 4
#define TS 32
#define RTS TS/WPT
__kernel void basic_multiply(const int M, const int N, const int K,
                      const __global double* A,
                      const __global double* B,
                      __global double* C) {
        
    const int row = get_local_id(0);
    const int col = get_local_id(1);
    const int i = TS*get_group_id(0) + row;
    const int j = TS*get_group_id(1) + col;
    
    __local double Asub[TS][TS];
    __local double Bsub[TS][TS];
 
    double acc[WPT];
    for (int w=0; w<WPT; w++) {
        acc[w] = 0.0;
    }
    
    const int numTiles = (K+TS-1)/TS;
    for (int t=0; t<numTiles; t++) {
        for (int w=0; w<WPT; w++) {
            const int tiled_i = TS*t + row;
            const int tiled_j = TS*t + col;
            if ( i < M && (tiled_j + w*RTS) < K ) {
              Asub[col + w*RTS][row] = A[(tiled_j + w*RTS)*M + i];
            }else{
              Asub[col + w*RTS][row] = 0.0;
            }
            if ( tiled_i < K && (j + w*RTS) < N ) {
              Bsub[col + w*RTS][row] = B[(j + w*RTS)*K + tiled_i];
            }else{
              Bsub[col + w*RTS][row] = 0.0;
            }
        }
        
        
        barrier(CLK_LOCAL_MEM_FENCE);
 
        for (int k=0; k<TS; k++) {
            for (int w=0; w<WPT; w++) {
                acc[w] += Asub[k][row] * Bsub[col + w*RTS][k];
            }
        }
 
        barrier(CLK_LOCAL_MEM_FENCE);
    }
 
    for (int w=0; w<WPT; w++) {
      if ( (j + w*RTS) < N && i < M ) {
        C[(j + w*RTS)*M + i] = acc[w];
      }
    }
}

#define WPT1 4
#define TS1 32
#define RTS1 TS1/WPT1
__kernel void multiply_self_transposed(const int M, const int N, const int K,
                      const __global double* A,
                      const __global double* B,
                      __global double* C) {
    
    const int row = get_local_id(0); 
    const int col = get_local_id(1); 
    const int i = TS1*get_group_id(0) + row; 
    const int j = TS1*get_group_id(1) + col; 
    const int jMin = TS1*get_group_id(1);
    const int iMax = TS1*get_group_id(0)  + get_local_size(0);
    
    __local double Asub[TS1][TS1];
    __local double Bsub[TS1][TS1];
 
    double acc[WPT1];
    for (int w=0; w<WPT1; w++) {
        acc[w] = 0.0;
    }
    
    
    const int numTiles = (K+TS1-1)/TS1;
    for (int t=0; t<numTiles; t++) {
        if(jMin <= iMax){
        
        for (int w=0; w<WPT1; w++) {
            const int tiled_i = TS1*t + row;
            const int tiled_j = TS1*t + col;
            if ( i < M && (tiled_j + w*RTS1) < K ) {
              Asub[col + w*RTS1][row] = A[(tiled_j + w*RTS1)*M + i];
            }else{
              Asub[col + w*RTS1][row] = 0.0;
            }
            if ( tiled_i < K && (j + w*RTS1) < N ) {
              Bsub[col + w*RTS1][row] = B[(j + w*RTS1)*K + tiled_i];
            }else{
              Bsub[col + w*RTS1][row] = 0.0;
            }
        }
        }
        
        barrier(CLK_LOCAL_MEM_FENCE);
 
        
        for (int k=0; k<TS1; k++) {
            for (int w=0; w<WPT1; w++) {
              if((j + w*RTS1)<=i){
                acc[w] += Asub[k][row] * Bsub[col + w*RTS1][k];
              }
            }
        }
 
        
        barrier(CLK_LOCAL_MEM_FENCE);
    }
 
    for (int w=0; w<WPT1; w++) {
      if ( (j + w*RTS1) < N && i < M  && (j + w*RTS1)<=i ) {
        C[(j + w*RTS1)*M + i] = acc[w];
        C[i*M + (j + w*RTS1)] = acc[w];
      }
    }
 }
)====="
