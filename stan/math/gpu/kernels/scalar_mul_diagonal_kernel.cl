R"(
#ifndef A
#define A(i, j)  A[j*rows+i]
#endif
__kernel void scalar_mul_diagonal(
      __global double *A,
      double scalar,
      unsigned int rows,
      unsigned int cols) {
  int i = get_global_id(0);
  if (i < rows && i < cols) {
   A(i, i) *= scalar;
  }
};
)"
