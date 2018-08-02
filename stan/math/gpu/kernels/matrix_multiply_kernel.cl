R"(
#ifndef WPT
#define WPT 8
#endif
#define RTS TS/WPT
/**
 * Matrix multiplication on the GPU
 *
 * @param M Number of rows for matrix A
 * @param N Number of cols for matrix A and rows for matrix B
 * @param K Number of cols for matrix B
 * @param[in] A the left matrix in matrix multiplication
 * @param[in] B the right matrix in matrix multiplication
 * @param[out] C the output matrix
 */
__kernel void matrix_multiply(const int M, const int N, const int K,
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
      
      Asub[col + w*RTS][row] = A[(tiled_j + w*RTS)*M + i];
      
      Bsub[col + w*RTS][row] = B[(j + w*RTS)*K + tiled_i];
      
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
   C[(j + w*RTS)*M + i] = acc[w];   
  }
  
}
)"
