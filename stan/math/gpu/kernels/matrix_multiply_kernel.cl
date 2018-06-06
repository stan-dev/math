R"(
#define WPT 4
#define RTS TS/WPT
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
};
)"
