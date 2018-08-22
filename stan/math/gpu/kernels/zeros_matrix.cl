#define STRINGIFY(src) #src

STRINGIFY(
/**
 * Stores zeros in the matrix on the GPU.
 * Supports writing zeroes to the lower and upper triangular or
 * the whole matrix.
 *
 * @param[out] A matrix
 * @param rows Number of rows for matrix A
 * @param cols Number of columns for matrix A
 * @param part optional parameter that describes where to assign zeros:
 *  LOWER - lower triangular
 *  UPPER - upper triangular
 * if the part parameter is not specified,
 * zeros are assigned to the whole matrix.
 *
 * @note  This kernel uses the helper macros available in helpers.cl.
 */
__kernel void zeros(__global write_only double *A, read_only unsigned int rows,
	read_only unsigned int cols, read_only unsigned int part) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < rows && j < cols) {
    if (part == LOWER && j < i) {
      A(i,j) = 0;
    } else if (part == UPPER && j > i) {
      A(i,j) = 0;
    } else if (part == ENTIRE) {
      A(i,j) = 0;
    }
  }
}
)
