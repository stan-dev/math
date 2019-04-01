#ifndef STAN_MATH_OPENCL_KERNELS_HELPERS_HPP
#define STAN_MATH_OPENCL_KERNELS_HELPERS_HPP
#ifdef STAN_OPENCL

#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

/*
 * Defines some helper macros for the kernels
 */
static const char* helpers =
    R"(
    // Matrix access helpers
  #ifndef A
  #define A(i,j) A[(j) * rows + (i)]
  #endif
  #ifndef B
  #define B(i,j) B[(j) * rows + (i)]
  #endif
  #ifndef C
  #define C(i,j) C[(j) * rows + (i)]
  #endif
    // Transpose
  #ifndef BT
  #define BT(i,j) B[(j) * cols + (i)]
  #endif
  #ifndef AT
  #define AT(i,j) A[(j) * cols + (i)]
  #endif
    // Moving between two buffers
  #ifndef src
  #define src(i,j) src[(j) * src_rows + (i)]
  #endif
  #ifndef dst
  #define dst(i,j) dst[(j) * dst_rows + (i)]
  #endif

  // The local memory column for each thread block
  #define THREAD_BLOCK_SIZE_COL THREAD_BLOCK_SIZE/WORK_PER_THREAD

  )";
}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
