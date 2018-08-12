#define STRINGIFY(src) #src

STRINGIFY(
/**
 * Copies a lower/upper triangular of a matrix to it's upper/lower.
 *
 * @param[in,out] A The matrix.
 * @param rows The number of rows in A.
 * @param cols The number of cols in A.
 * @param copy_direction A value of zero or one specifying
 *  which direction to copy
 *  LOWER_TO_UPPER: 1
 *  UPPER_TO_LOWER: 0
 *
 * @note Used in mat/gpu/copy_triangular_transposed.hpp.
 *  This kernel uses the helper macros available in helpers.cl.
 */
__kernel void copy_triangular_transposed(__global read_write double *A,
	read_only unsigned int rows, read_only unsigned int cols,
	read_only unsigned int copy_direction) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < rows && j < cols) {
    if (copy_direction == LOWER_TO_UPPER && i > j) {
      A(j, i) = A(i, j);
    } else if (copy_direction == UPPER_TO_LOWER && i > j) {
      A(i, j) = A(j, i);
    }
  }
}
);
