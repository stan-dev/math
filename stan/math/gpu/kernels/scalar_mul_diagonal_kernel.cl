R"(
#ifndef A
#define A(i, j)  A[j*rows+i]
#endif
/**
 * Matrix multiplication of the form A*A^T on the GPU
 *
 * @param[in, out] A matrix A
 * @param scalar the value with which to multiply the diagonal of A
 * @param rows the number of rows in A
 * @param cols the number of columns in A
 */
__kernel void scalar_mul_diagonal(
      __global double *A,
      double scalar,
      unsigned int rows,
      unsigned int cols) {
  int i = get_global_id(0);
  if (i < rows && i < cols) {
   A(i, i) *= scalar;
  }
};
)"
