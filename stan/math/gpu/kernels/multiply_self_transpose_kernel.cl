R"(
#ifndef WORK_PER_THREAD_MULT_SELF_TRANS
#define WORK_PER_THREAD_MULT_SELF_TRANS 4
#endif
#define WORKGROUP_SIZE_MULT_SELF_TRANS_COL WORKGROUP_SIZE_MULT_SELF_TRANS/WORK_PER_THREAD_MULT_SELF_TRANS
#ifndef B
#define B(i, j)  B[j * M + i]
#endif

/**
 * Matrix multiplication of the form A*A^T on the GPU
 *
 * @param[in] A matrix A
 * @param[out] B the output matrix
 * @param M Number of rows for matrix A
 * @param N Number of cols for matrix A and the number of rows for matrix A^T 
 */
__kernel void multiply_self_transpose(
            const __global double* A,
            __global double* B,
            const int M,
            const int N) {
  // thread index inside the workgroup
  const int workgroup_row = get_local_id(0); 
  const int workgroup_col = get_local_id(1); 
  
  // global thread index
  const int i = WORKGROUP_SIZE_MULT_SELF_TRANS*get_group_id(0) + workgroup_row; 
  const int j = WORKGROUP_SIZE_MULT_SELF_TRANS*get_group_id(1) + workgroup_col; 
  
  // indexes that determine the last indexes that need to compute
  // in order to remove the unnecesary multiplications in the special multiplication of A*A^T
  const int jMin = WORKGROUP_SIZE_MULT_SELF_TRANS*get_group_id(1);
  const int iMax = WORKGROUP_SIZE_MULT_SELF_TRANS*get_group_id(0)  + get_local_size(0);
  
  // local memory
  __local double A_local[WORKGROUP_SIZE_MULT_SELF_TRANS][WORKGROUP_SIZE_MULT_SELF_TRANS];
  __local double B_local[WORKGROUP_SIZE_MULT_SELF_TRANS][WORKGROUP_SIZE_MULT_SELF_TRANS];

  double acc[WORK_PER_THREAD_MULT_SELF_TRANS];
  for (int w = 0; w < WORK_PER_THREAD_MULT_SELF_TRANS; w++) {
    acc[w] = 0.0;
  }
  
  const int numTiles = (N+WORKGROUP_SIZE_MULT_SELF_TRANS-1)/WORKGROUP_SIZE_MULT_SELF_TRANS;
  //iterate over all tiles
  for (int t = 0; t < numTiles; t++) {
    // in each tile
    const int tiled_i = WORKGROUP_SIZE_MULT_SELF_TRANS*t + workgroup_row;
    const int tiled_j = WORKGROUP_SIZE_MULT_SELF_TRANS*t + workgroup_col;
    // if the data needs to be loaded to local memory
    if(jMin <= iMax){
      // each thread copies WORK_PER_THREAD_MULT_SELF_TRANS values to the local memory
      for (int w = 0; w < WORK_PER_THREAD_MULT_SELF_TRANS; w++) {
        A_local[workgroup_col + w*WORKGROUP_SIZE_MULT_SELF_TRANS_COL][workgroup_row] = A[i + (tiled_j + w*WORKGROUP_SIZE_MULT_SELF_TRANS_COL)*M];
        B_local[workgroup_col + w*WORKGROUP_SIZE_MULT_SELF_TRANS_COL][workgroup_row] = A[(j + w*WORKGROUP_SIZE_MULT_SELF_TRANS_COL) + tiled_i*M];      
      }
    }
    // wait till all tile values are loaded to the local memory
    barrier(CLK_LOCAL_MEM_FENCE);
    // multiply the tile producs
    for (int k = 0; k < WORKGROUP_SIZE_MULT_SELF_TRANS; k++) {
      // each thread multiplies WORK_PER_THREAD_MULT_SELF_TRANS values
      for (int w = 0;w < WORK_PER_THREAD_MULT_SELF_TRANS; w++) {
        if(jMin <= iMax){
          if((j + w*WORKGROUP_SIZE_MULT_SELF_TRANS_COL)<=i){
            acc[w] += A_local[k][workgroup_row] * B_local[workgroup_col + w*WORKGROUP_SIZE_MULT_SELF_TRANS_COL][k];
          }
        }
      }
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  // save the values
  for (int w=0; w<WORK_PER_THREAD_MULT_SELF_TRANS; w++) {
    // each thread saves WORK_PER_THREAD_MULT_SELF_TRANS values
    if ( (j + w*WORKGROUP_SIZE_MULT_SELF_TRANS_COL) <= i ) {
      B[i + (j + w*WORKGROUP_SIZE_MULT_SELF_TRANS_COL)*M] = acc[w];
      B[(j + w*WORKGROUP_SIZE_MULT_SELF_TRANS_COL) + i*M] = acc[w];
    }
  }
};
)"
