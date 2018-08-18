#define STRINGIFY(src) #src

STRINGIFY(
/**
 * Makes an identity matrix on the GPU
 *
 * @param[in,out] A The identity matrix output.
 * @param rows The number of rows for A.
 * @param cols The number of cols for A.
 *
 * @note Used in math/gpu/identity_opencl.hpp.
 *  This kernel uses the helper macros available in helpers.cl.
 */
__kernel void identity(__global write_only double *A,
	read_only unsigned int rows, read_only unsigned int cols) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < rows && j < cols) {
    if (i == j) {
      A(i, j) = 1.0;
    } else {
      A(i, j) = 0.0;
    }
  }
}
)
