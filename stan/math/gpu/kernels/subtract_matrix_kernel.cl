R"(
#ifndef A
#define A(i, j)  A[j * rows + i]
#endif
#ifndef B
#define B(i, j)  B[j * rows + i]
#endif
#ifndef C
#define C(i, j)  C[j * rows + i]
#endif

/**
 * Matrix subtraction on the GPU Subtracts the second matrix
 * from the first matrix and stores the result in the third matrix (C=A-B).
 *
 * @param C (write) The output matrix.
 * @param B (read) RHS input matrix.
 * @param A (read) LHS input matrix.
 * @param rows The number of rows for matrix A.
 * @param cols The number of columns for matrix A.
 *
 * @note Used in math/gpu/subtract_opencl.hpp
 */
__kernel void subtract(__global double *C, __global double *A,
  __global double *B, unsigned int rows, unsigned int cols) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < rows && j < cols) {
    C(i, j) = A(i, j) - B(i, j);
  }
};)"
