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

	// Options for lower or upper triangular as well as full
	#ifndef LOWER
	#define LOWER 0
	#endif
	#ifndef UPPER
	#define UPPER 1
	#endif
	#ifndef ENTIRE
	#define ENTIRE 2
	#endif

	// Options for flipping lower to upper and vice versa
	#ifndef UPPER_TO_LOWER
	#define UPPER_TO_LOWER 0
	#endif
	#ifndef LOWER_TO_UPPER
	#define LOWER_TO_UPPER 1
	#endif

	// Thread block sizes
  #ifndef WORK_PER_THREAD_MULT
  #define WORK_PER_THREAD_MULT 8
  #endif
  #ifndef THREAD_BLOCK_SIZE
  #define THREAD_BLOCK_SIZE 32
  #endif
  #ifndef THREAD_BLOCK_SIZE_MULT_COL
  #define THREAD_BLOCK_SIZE_MULT_COL THREAD_BLOCK_SIZE/WORK_PER_THREAD_MULT
  #endif
  #ifndef THREAD_BLOCK_SIZE
  #define THREAD_BLOCK_SIZE 32
  #endif
  #ifndef WORK_PER_THREAD_MULT_SELF_TRANS
  #define WORK_PER_THREAD_MULT_SELF_TRANS 4
  #endif
  #define THREAD_BLOCK_SIZE_MULT_SELF_TRANS_COL \
    THREAD_BLOCK_SIZE / WORK_PER_THREAD_MULT_SELF_TRANS
  )";
}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
