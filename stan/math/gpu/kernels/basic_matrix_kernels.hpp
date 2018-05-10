#ifndef STAN_MATH_GPU_BASIC_MATRIX_KERNELS_HPP
#define STAN_MATH_GPU_BASIC_MATRIX_KERNELS_HPP
#ifdef STAN_OPENCL
#include <string>

/**
 *  @file stan/math/gpu/kernels/basic_matrix_kernels.hpp
 *  @brief Kernel sources for basic matrix operations on the gpu:
 *    copy, copy lower/upper triangular, copy triangular transposed,
 *    copy submatrix, init to zeros, init to identity,
 *    add, subtract, transpose
 */

namespace stan {
namespace math {

std::string copy_matrix_kernel
    = "#define A(i,j)  A[j*rows+i] \n"
      "#define B(i,j)  B[j*rows+i] \n"
      " __kernel void copy( \n"
      "      __global double *A, \n"
      "      __global double *B, \n"
      "      unsigned int rows, \n"
      "      unsigned int cols) { \n"
      "  int i = get_global_id(0); \n"
      "  int j = get_global_id(1); \n"
      "  if ( i < rows && j < cols ) { \n"
      "   B(i,j) = A(i,j); \n"
      "  }\n"
      "}\n";
}
}  // namespace stan
#endif
#endif
