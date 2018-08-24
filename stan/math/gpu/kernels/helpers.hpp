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
  // Helper macros
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
  #ifndef WORK_PER_THREAD_MULT
  #define WORK_PER_THREAD_MULT 8
  #endif
  #ifndef WORKGROUP_SIZE_MULT
  #define WORKGROUP_SIZE_MULT 32
  #endif
  #ifndef WORKGROUP_SIZE_MULT_COL
    #define WORKGROUP_SIZE_MULT_COL WORKGROUP_SIZE_MULT/WORK_PER_THREAD_MULT
  #endif
  #ifndef WORKGROUP_SIZE_MULT_SELF_TRANS
  #define WORKGROUP_SIZE_MULT_SELF_TRANS 32
  #endif
  #ifndef WORK_PER_THREAD_MULT_SELF_TRANS
    #define WORK_PER_THREAD_MULT_SELF_TRANS 4
  #endif
  #define WORKGROUP_SIZE_MULT_SELF_TRANS_COL \
    WORKGROUP_SIZE_MULT_SELF_TRANS / WORK_PER_THREAD_MULT_SELF_TRANS
  )";
}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
