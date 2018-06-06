R"(
#ifndef A
#define A(i, j)  A[j*rows+i]
#endif
#ifndef B
#define B(i, j)  B[j*rows+i]
#endif
__kernel void scalar_mul(
      __global double *A,
      __global double *B,
      double scalar,
      unsigned int rows,
      unsigned int cols) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < rows && j < cols) {
   A(i, j) = B(i, j)*scalar;
  }
};)"
