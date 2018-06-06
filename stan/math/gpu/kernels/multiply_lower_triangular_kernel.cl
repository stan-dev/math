R"(
#define WPT1 4
#define RTS1 TS1/WPT1
__kernel void multiply_lower_triangular(
            const int M,
            const int N,
            const int K,
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
        if((k+t*TS1)>=i)
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
};  
)"
