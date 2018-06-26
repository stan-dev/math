R"(
#ifndef A
#define A(i, j)  A[j * rows + i]
#endif

/**
 * Makes an identity matrix on the GPU
 *
 * @param A (write) The identity matrix output.
 * @param rows The number of rows for A.
 * @param cols The number of cols for A.
 *
 * @note Used in math/gpu/identity_opencl.hpp
 */
__kernel void identity(__global double *A, unsigned int rows,
  unsigned int cols) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < rows && j < cols) {
    if (i == j) {
      A(i, j) = 1.0;
    } else {
      A(i, j) = 0.0;
    }
  }
};)"
