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
#define A(i, j)  A[j * rows + i]
#endif
/**
 * Stores zeros in the matrix on the GPU.
 * Supports writing zeroes to the lower and upper triangular or
 * the whole matrix.
 *
 * @param A (write) matrix
 * @param rows Number of rows for matrix A
 * @param cols Number of columns for matrix A
 * @param part optional parameter that describes where to assign zeros:
 *  LOWER - lower triangular
 *  UPPER - upper triangular
 * if the part parameter is not specified,
 * zeros are assigned to the whole matrix.
 *
 */
__kernel void zeros(__global double *A, unsigned int rows, unsigned int cols,
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
