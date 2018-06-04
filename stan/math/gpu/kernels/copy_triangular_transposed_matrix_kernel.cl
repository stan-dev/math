R"(
#ifndef LOWER_TO_UPPER
#define LOWER_TO_UPPER 1
#endif

#ifndef UPPER_TO_LOWER
#define UPPER_TO_LOWER 0
#endif

#ifndef A
#define A(i, j)  A[j*rows+i]
#endif

#ifndef AT
#define AT(i, j)  A[j*cols+i] 
#endif
__kernel void copy_triangular_transposed(
      __global double *A,
      unsigned int rows,
      unsigned int cols,
      unsigned int copy_direction) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < rows && j < cols) {
    if (copy_direction==LOWER_TO_UPPER && j > i) {
      AT(j, i) = A(i, j);
    } else if (copy_direction==UPPER_TO_LOWER && j > i) {
      A(i, j) = AT(j, i);
    }
  }
};)"
