#ifndef STAN_MATH_GPU_KERNELS_HELPERS_HPP
#define STAN_MATH_GPU_KERNELS_HELPERS_HPP
#ifdef STAN_OPENCL

#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

/*
 * Defines some helper macros for the kernels
 */
std::string helpers =
    R"(
  // Matrix access helpers
  #ifndef A
  #define A(i,j) A[j * rows + i]
  #endif
  #ifndef B
  #define B(i,j) B[j * rows + i]
  #endif
  #ifndef C
  #define C(i,j) C[j * rows + i]
  #endif
  #ifndef BT
  #define BT(i,j) B[j * cols + i]
  #endif
  #ifndef src
  #define src(i,j) src[j * src_rows + i]
  #endif
  #ifndef dst
  #define dst(i,j) dst[j * dst_rows + i]
  #endif

	// Thread block sizes
  #ifndef WORK_PER_THREAD
  #define WORK_PER_THREAD 8
  #endif
  #ifndef THREAD_BLOCK_SIZE
  #define THREAD_BLOCK_SIZE 32
  #endif
  #ifndef THREAD_BLOCK_SIZE_COL
  #define THREAD_BLOCK_SIZE_COL THREAD_BLOCK_SIZE/WORK_PER_THREAD
  #endif
  )";
}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
