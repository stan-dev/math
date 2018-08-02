R"(
#ifndef A
#define A(i, j)  A[j*M+i]
#endif
/**
 * Matrix multiplication of the form A*A^T on the GPU
 *
 * @param[in, out] A matrix A
 * @param scalar the value with which to multiply the diagonal of A
 * @param M the number of rows in A
 * @param N the number of columns in A
 */
__kernel void scalar_mul_diagonal(
      __global double *A,
      const double scalar,
      const unsigned int M,
      const unsigned int N) {
  int i = get_global_id(0);
  if (i < M && i < N) {
   A(i, i) *= scalar;
  }
};
)"
