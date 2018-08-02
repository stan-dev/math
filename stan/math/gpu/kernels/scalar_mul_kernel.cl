R"(
#ifndef A
#define A(i, j)  A[j*M+i]
#endif
#ifndef B
#define B(i, j)  B[j*M+i]
#endif
__kernel void scalar_mul(
      __global double *A,
      const __global double *B,
      const double scalar,
      const unsigned int M,
      const unsigned int N) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < M && j < N) {
   A(i, j) = B(i, j)*scalar;
  }
};)"
