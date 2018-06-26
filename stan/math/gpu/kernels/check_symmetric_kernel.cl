R"(
#ifndef A
#define A(i, j)  A[j * rows + i]
#endif
/**
 * Check if the <code>matrix_gpu</code> is symmetric
 *
 * @param A The matrix to check.
 * @param rows The number of rows in matrix A.
 * @param cols The number of columns in matrix A.
 * @param flag (write) the flag to be written to if any diagonal is zero.
 *
 * @note Kernel for stan/math/gpu/err/check_symmetric.hpp
 */
__kernel void is_symmetric(__global double *A, int rows, int cols,
  __global int *flag, double tolerance) {
  const int i = get_global_id(0);
  const int j = get_global_id(1);
  if (i < rows && j < cols) {
    double diff = fabs(A(i, j) - A(j, i));
    if (diff > tolerance) {
      flag[0] = 0;
    }
  }
};)"
