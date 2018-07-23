R"(
#ifndef LOWER_TO_UPPER
#define LOWER_TO_UPPER 1
#endif

#ifndef UPPER_TO_LOWER
#define UPPER_TO_LOWER 0
#endif

#ifndef A
#define A(i, j)  A[j * rows + i]
#endif

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
 * @note Used in mat/gpu/copy_triangular_transposed.hpp
 */
__kernel void copy_triangular_transposed(__global double *A, unsigned int rows,
  unsigned int cols, unsigned int copy_direction) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < rows && j < cols) {
    if (copy_direction == LOWER_TO_UPPER && i > j) {
      A(j, i) = A(i, j);
    } else if (copy_direction == UPPER_TO_LOWER && i > j) {
      A(i, j) = A(j, i);
    }
  }
};)"
