R"(
#ifndef A
#define A(i, j)  A[j * rows + i]
#endif
/**
 * Check if the <code>matrix_gpu</code> has zeros on the diagonal
 *
 * @param A Matrix to check.
 * @param rows The number of rows for A.
 * @param cols The number of cols of A.
 * @param flag (write) the flag to be written to if any diagonal is zero.
 *
 * @note Kernel for stan/math/gpu/err/check_diagonal_zeros.hpp
 */
__kernel void is_zero_on_diagonal(__global double *A, int rows, int cols,
  __global int *flag) {
  const int i = get_global_id(0);
  if (i < rows && i < cols) {
    if (A(i, i) == 0) {
      flag[0] = 1;
    }
  }
};)"
