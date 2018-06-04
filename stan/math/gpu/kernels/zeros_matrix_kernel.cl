R"(
#ifndef LOWER
#define LOWER 0
#endif

#ifndef UPPER
#define UPPER 1
#endif

#ifndef ALL
#define ALL 2
#endif

#ifndef A
#define A(i, j)  A[j*rows+i]
#endif

#ifndef BT
#define BT(i, j)  B[j*cols+i] 
#endif

__kernel void zeros(
        __global double *A,
        unsigned int rows,
        unsigned int cols,
        unsigned int part) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < rows && j < cols) {
    if (part == LOWER && j < i) {
      A(i,j) = 0;
    } else if (part == UPPER && j > i) {
      A(i,j) = 0;
    } else if (part == ALL) {
      A(i,j) = 0;
    }
  }
};)"
