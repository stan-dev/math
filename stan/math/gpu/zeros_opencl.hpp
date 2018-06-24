#ifndef STAN_MATH_GPU_ZEROS_OPENCL_HPP
#define STAN_MATH_GPU_ZEROS_OPENCL_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/matrix_gpu.hpp>
#include <CL/cl.hpp>

/** @file stan/math/gpu/basic_matrix_gpu.hpp
 *    @brief basic_matrix_gpu - basic matrix operations:
 *    copy, copy lower/upper triangular, copy triangular transposed,
 *    copy submatrix, init matrix with zeros, create an identity matrix,
 *    add, subtract, transpose
 */

namespace stan {
namespace math {
enum triangularity { LOWER = 0, UPPER = 1, NONE = 2 };
/**
 * Stores zeros in the matrix on the GPU.
 * It support writing the zeros to lower
 * triangular, upper triangular or the
 * whole matrix.
 *
 * @param A matrix
 * @param part optional parameter
 * that describes where to assign zeros:
 * LOWER - lower triangular
 * UPPER - upper triangular
 * if the parameter is not specified,
 * zeros are assigned to the whole matrix.
 *
 */
inline void zeros(matrix_gpu& A, triangularity part = NONE) {
  if (A.size() == 0)
    return;
  cl::Kernel kernel = opencl_context.get_kernel("zeros");
  cl::CommandQueue cmdQueue = opencl_context.queue();
  try {
    kernel.setArg(0, A.buffer());
    kernel.setArg(1, A.rows());
    kernel.setArg(2, A.cols());
    kernel.setArg(3, part);
    cmdQueue.enqueueNDRangeKernel(kernel, cl::NullRange,
                                  cl::NDRange(A.rows(), A.cols()),
                                  cl::NullRange, NULL, NULL);
  } catch (const cl::Error& e) {
    check_opencl_error("zeros", e);
  }
}
}  // namespace math
}  // namespace stan

#endif
#endif
