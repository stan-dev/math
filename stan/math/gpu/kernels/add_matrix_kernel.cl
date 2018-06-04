R"(
#ifndef A
#define A(i, j)  A[j*rows+i]
#endif
#ifndef B
#define B(i, j)  C[j*rows+i]
#endif
#ifndef C
#define C(i, j)  C[j*rows+i]
#endif
__kernel void add(
      __global double *C,
      __global double *A,
      __global double *B,
      unsigned int rows,
      unsigned int cols) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < rows && j < cols) {
    C(i, j) = A(i, j) + B(i, j);
  }
};)"
