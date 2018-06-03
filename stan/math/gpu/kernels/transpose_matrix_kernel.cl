R"(
#ifndef A
#define A(i, j)  A[j*rows+i]
#endif
#ifndef BT
#define BT(i, j)  B[j*cols+i] 
#endif

__kernel void transpose(
        __global double *B,
        __global double *A,
        int rows,
        int cols ) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < rows && j < cols) {
    BT(j, i) = A(i, j);
  }
};)"
