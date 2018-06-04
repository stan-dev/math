R"(
#ifndef A
#define A(i, j)  A[j*rows+i]
#endif
#ifndef BT
#define BT(i, j)  B[j*cols+i] 
#endif

__kernel void identity(
      __global double *A,
      unsigned int rows,
      unsigned int cols) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < rows && j < cols) {
    if (i == j) {
      A(i, j) = 1.0;
    } else {
      A(i, j) = 0.0;
    }
  }
};)"
