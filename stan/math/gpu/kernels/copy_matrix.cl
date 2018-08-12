#define STRINGIFY(src) #src

STRINGIFY(
/**
 * Copy one matrix to another
 * @param[in] A The matrix to copy.
 * @param[out] B The matrix to copy A to.
 * @param rows The number of rows in A.
 * @param cols The number of cols in A.
 *
 * @note Kernel used in math/gpu/matrix_gpu.hpp
 */
__kernel void copy(__global read_only double *A, __global write_only double *B,
	read_only unsigned int rows, read_only unsigned int cols) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < rows && j < cols) {
    B(i, j) = A(i, j);
  }
}
);
