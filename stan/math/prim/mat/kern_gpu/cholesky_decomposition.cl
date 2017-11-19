R"=====(
#ifdef cl_khr_fp64
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
    #error "Double not supported by OpenCL implementation."
#endif

__kernel void cholesky_block(
                __global double *V,
                __global double *d,
                int n) {
  int f = get_local_id(0);
  for (int j = 0; j <n; j++) { 
      double ss = 0;
      for (int k = 0; k < j; k++) {
          ss += V[j * n + k] * V[j * n + k];
      }
      V[j * n + j] = sqrt(d[j * n + j] - ss);
      barrier(CLK_LOCAL_MEM_FENCE);
      int i=f;
      if ( i >= (j+1) && i < n ){
          double s = 0;
          for (int k = 0; k < j; k++) {
              s += V[i * n + k] * V[j * n + k];
          }
          V[i * n + j] = (1.0 / V[j * n + j] * (d[i * n + j] - s));
      }
  }
}

#define TS3 16
__kernel void cholesky_left_update(
                __global double* l,
                __global double* ap,
                __global double* temp,
                int offset,
                int block,
                int M,
                int max_threads) {
  int k;
  const int row = get_local_id(0);
  const int col = get_local_id(1);
  const int i = TS3*get_group_id(0) + row;
  const int j = TS3*get_group_id(1) + col;
  __local double Asub[TS3][TS3];
  __local double Bsub[TS3][TS3];
  double sum = 0;
  const int numTiles = (block+TS3-1)/TS3;
  for (int t=0; t < numTiles; t++) {
    const int tiledRow = TS3*t + row;
    const int tiledCol = TS3*t + col;
    if ( i < max_threads && tiledCol < block &&
        (i+offset+block) < M && (offset+tiledCol) < M) {
      Asub[col][row] = ap[(i+offset+block)*M+offset+tiledCol];
    } else {
      Asub[col][row] = 0.0;
    }
    if ( j < block && tiledRow < block ) {
      Bsub[row][col] = temp[j*block+tiledRow];
    } else {
      Bsub[row][col] = 0.0;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    for (k = 0; k < TS3; k++) {
      sum += Asub[k][row]*Bsub[k][col];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  if ( i < max_threads && j < block && i < M && j < M ) {
    l[i*block+j] = sum;
  }  
}

__kernel void cholesky_zero(
                __global double* A,
                int M) {
  const int i = get_global_id(0);
  const int j = get_global_id(1);
  if ( j > i ) {
    A[i*M+j] = 0.0;
  }
}
)====="
