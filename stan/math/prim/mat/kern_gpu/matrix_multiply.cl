R"=====(
#ifdef cl_khr_fp64
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
    #error "Double not supported by OpenCL implementation."
#endif

#define a(i,j)  a[i*cols+j]
#define b(i,j)  b[i*cols+j]
#define c(i,j)  c[i*cols+j]

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

#define A(i,j)  A[i*K+j]
#define B(i,j)  B[i*N+j]
#define C(i,j)  C[i*N+j]

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
        
        const int tiledRow = TS*t + local_i;
        const int tiledCol = TS*t + local_j;

        if ( tiledRow < K && j < N ) {
         Asub[local_j][local_i] = B(tiledRow,j);
        } else {
         Asub[local_j][local_i] = 0.0;
        }

        if ( tiledCol < K && i < M ) {
         Bsub[local_j][local_i] = A(i,tiledCol);
        } else {
         Bsub[local_j][local_i] = 0.0;
        }

        barrier(CLK_LOCAL_MEM_FENCE);

        for (int k=0; k < TS; k++) {
            acc += Bsub[k][local_i] * Asub[local_j][k];
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
    // Thread identifiers
    const int local_i = get_local_id(0);  // Local row ID (max: TS)
    const int local_j = get_local_id(1);  // Local col ID (max: TS)
    const int i = TS11*get_group_id(0) + local_i;  // Row ID of C (0..M)
    const int j = TS11*get_group_id(1) + local_j;  // Col ID of C (0..N)
    const int jMin = TS11*get_group_id(1);
    const int iMax = TS11*get_group_id(0)  + get_local_size(0);

    // Local memory to fit a tile of TS*TS elements of A and B
    __local double Asub[TS11][TS11];
    __local double Bsub[TS11][TS11];

    // Initialise the accumulation register
    double acc = 0.0;

    // Loop over all tiles
    const int numTiles = (K+TS11-1)/TS11;
    for (int t=0; t < numTiles; t++) {
        // Load one tile of A and B into local memory
        const int tiledRow = TS11*t + local_i;
        const int tiledCol = TS11*t + local_j;
        if(jMin <= iMax){
          if ( tiledRow < K && j < N ) {
           Asub[local_j][local_i] = B(tiledRow,j);
          } else {
           Asub[local_j][local_i] = 0.0;
          }
          if ( tiledCol < K && i < M ) {
           Bsub[local_j][local_i] = A(i,tiledCol);
          } else {
           Bsub[local_j][local_i] = 0.0;
          }
        }
        // Synchronise to make sure the tile is loaded
        barrier(CLK_LOCAL_MEM_FENCE);

        // Perform the computation for a single tile
        if( j <= i ){
          for (int k=0; k < TS11; k++) {
              acc += Bsub[k][local_i] * Asub[local_j][k];
          }
        }

        // Synchronise before loading the next tile
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // Store the final result in C
    if ( j < N && i < M && j <= i ) {
     C(i,j) = acc;
     C(j,i) = acc;
    }
 }
)====="
