#define STRINGIFY(src) #src

STRINGIFY(
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
__kernel void add(__global write_only double *C, __global read_only double *A,
	__global read_only double *B, read_only unsigned int rows,
	 read_only unsigned int cols) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < rows && j < cols) {
    C(i, j) = A(i, j) + B(i, j);
  }
};
);
