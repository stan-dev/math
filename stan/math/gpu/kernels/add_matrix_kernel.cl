#define STRINGIFY(src) #src

STRINGIFY(
//#ifndef A
//#define A(i, j)  A[j * rows + i]
//#endif
//#ifndef B
//#define B(i, j)  B[j * rows + i]
//#endif
//#ifndef C
//#define C(i, j)  C[j * rows + i]
//#endif
/**
 * Matrix addition on the GPU
 *
 * @param[out] C Output matrix.
 * @param[in] A LHS of matrix addition.
 * @param[in] B RHS of matrix addition.
 * @param rows Number of rows for matrix A.
 * @param cols Number of rows for matrix B.
 *
 */
__kernel void add(__global double *C, __global double *A, __global double *B,
  unsigned int rows, unsigned int cols) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < rows && j < cols) {
    C[j * rows + i] = A[j * rows + i] + B[j * rows + i];
  }
}
);
