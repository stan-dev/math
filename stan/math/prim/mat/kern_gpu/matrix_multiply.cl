R"=====(
#ifdef cl_khr_fp64
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
    #error "Double precision floating point not supported by OpenCL implementation."
#endif

__kernel void scalar_mul_diagonal(
          __global double *a,
          double scalar,
          unsigned int M,
          unsigned int N) 
{ 
    int i = get_global_id(0); 
    if (i < M && i < N){ 
     a[i*N+i] *= scalar  ;
    } 
};
      
      
__kernel void scalar_mul(
          __global double *a,
          __global double *b,
          double scalar,
          unsigned int M,
          unsigned int N) 
{ 
    int i = get_global_id(0); 
    int j = get_global_id(1); 
    if (i < M && j < N){ 
     a[i*N+j] = b[i*N+j]*scalar; 
    } 
};

      
 #define TS 16 
__kernel void basic_multiply(const int M, const int N, const int K, 
                      const __global double* A, 
                      const __global double* B, 
                      __global double* C) { 
      
    // Thread identifiers 
    const int row = get_local_id(0); // Local row ID (max: TS)  
    const int col = get_local_id(1); // Local col ID (max: TS)  
    const int globalRow = TS*get_group_id(0) + row; // Row ID of C (0..M) 
    const int globalCol = TS*get_group_id(1) + col; // Col ID of C (0..N) 
  
    // Local memory to fit a tile of TS*TS elements of A and B  
    __local double Asub[TS][TS];  
    __local double Bsub[TS][TS];  
  
    // Initialise the accumulation register 
    double acc = 0.0; 
    
    // Loop over all tiles 
    const int numTiles = (K+TS-1)/TS;
    for (int t=0; t<numTiles; t++) {
 
        // Load one tile of A and B into local memory 
        const int tiledRow = TS*t + row; 
        const int tiledCol = TS*t + col; 
        if (tiledRow<K && globalCol<N){  
         Asub[col][row] = B[tiledRow*N + globalCol];  
        }else{
         Asub[col][row] = 0.0;  
        }
        if (tiledCol<K && globalRow<M){  
         Bsub[col][row] = A[globalRow*K + tiledCol];  
        }else{
         Bsub[col][row] = 0.0;  
        }
 
        // Synchronise to make sure the tile is loaded 
        barrier(CLK_LOCAL_MEM_FENCE);  
 
        // Perform the computation for a single tile 
        for (int k=0; k<TS; k++) { 
            acc += Bsub[k][row] * Asub[col][k];  
        }  
 
        // Synchronise before loading the next tile  
        barrier(CLK_LOCAL_MEM_FENCE);  
    }  
 
    // Store the final result in C 
    if (globalCol<N && globalRow<M){  
     C[globalRow*N + globalCol] = acc;  
    } 
 } 
)====="
