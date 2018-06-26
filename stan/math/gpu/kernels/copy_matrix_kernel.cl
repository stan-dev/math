R"(
#ifndef A
#define A(i, j)  A[j * rows + i]
#endif

#ifndef B
#define B(i, j)  B[j * rows + i]
#endif
/**
 * Copy one matrix to another
 * @param A (read) The matrix to copy.
 * @param B (write) The matrix to copy A to.
 * @param rows The number of rows in A.
 * @param cols The number of cols in A.
 *
 * @note Kernel used in math/gpu/matrix_gpu.hpp
 */
__kernel void copy(__global double *A, __global double *B, unsigned int rows,
  unsigned int cols) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < rows && j < cols) {
    B(i, j) = A(i, j);
  }
};)"
