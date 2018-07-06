R"(
#ifndef A
#define A(i, j)  A[j * rows + i]
#endif
#ifndef BT
#define BT(i, j)  B[j * cols + i]
#endif
/**
 * Takes the transpose of the matrix on the GPU.
 *
 * @param B (write) The output matrix to hold transpose of A.
 * @param A (read) The input matrix to transpose into B.
 * @param rows The number of rows for A.
 * @param cols The number of columns for A.
 *
 */
__kernel void transpose(__global double *B, __global double *A, int rows,
  int cols ) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < rows && j < cols) {
    BT(j, i) = A(i, j);
  }
};)"
