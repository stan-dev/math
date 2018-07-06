R"(
#ifndef A
#define A(i, j)  A[j * rows + i]
#endif

#ifndef B
#define B(i, j)  B[j * rows + i]
#endif
/**
 * Copies the lower or upper
 * triangular of the source matrix to
 * the destination matrix.
 * Both matrices are stored on the GPU.
 *
 * @param A (write) Output matrix to copy triangular to.
 * @param B (read) The matrix to copy the triangular from.
 * @param rows The number of rows of B.
 * @param cols The number of cols of B.
 * @param lower_upper enum to describe
 * which part of the matrix to copy:
 *  LOWER - copies the lower triangular
 *  UPPER - copes the upper triangular
 *
 * @note Used in math/gpu/copy_triangular_opencl.hpp
 */
__kernel void copy_triangular(__global double *A, __global double *B,
  unsigned int rows, unsigned int cols, unsigned int lower_upper) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < rows && j < cols) {
    if (!lower_upper && j <= i) {
      A(i, j) = B(i, j);
    } else if (!lower_upper) {
      A(i, j) = 0;
    } else if (lower_upper && j >= i) {
      A(i, j) = B(i, j);
    } else if (lower_upper && j < i) {
      A(i, j) = 0;
    }
  }
};)"
