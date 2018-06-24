#ifndef STAN_MATH_GPU_COPY_SUBMATRIX_HPP
#define STAN_MATH_GPU_COPY_SUBMATRIX_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/matrix_gpu.hpp>
#include <CL/cl.hpp>

namespace stan {
namespace math {
/**
 * Copies a submatrix of the source matrix to
 * the destination matrix. The submatrix to copy
 * starts at (src_offset_rows, src_offset_cols)
 * and is of size size_rows x size_cols.
 * The submatrix is copied to the
 * destination matrix starting at
 * (dst_offset_rows, dst_offset_cols)
 *
 * @param src the source matrix
 * @param dst the destination submatrix
 * @param src_offset_rows the offset row in src
 * @param src_offset_cols the offset column in src
 * @param dst_offset_rows the offset row in dst
 * @param dst_offset_cols the offset column in dst
 * @param size_rows the number of rows in the submatrix
 * @param size_cols the number of columns in the submatrix
 *
 * @throw <code>std::invalid_argument</code> if
 *
 *
 */
inline void copy_submatrix(matrix_gpu& src, matrix_gpu& dst,
                           int src_offset_rows, int src_offset_cols,
                           int dst_offset_rows, int dst_offset_cols,
                           int size_rows, int size_cols) {
  if (size_rows == 0 || size_cols == 0)
    return;
  if ((src_offset_rows + size_rows) > src.rows()
      || (src_offset_cols + size_cols) > src.cols()) {
    domain_error("copy_submatrix", "submatrix in src", " is out of bounds", "");
  }
  if ((dst_offset_rows + size_rows) > dst.rows()
      || (dst_offset_cols + size_cols) > dst.cols()) {
    domain_error("copy_submatrix", "submatrix in src", " is out of bounds", "");
  }
  cl::Kernel kernel = opencl_context.get_kernel("copy_submatrix");
  cl::CommandQueue cmdQueue = opencl_context.queue();
  try {
    kernel.setArg(0, src.buffer());
    kernel.setArg(1, dst.buffer());

    kernel.setArg(2, src_offset_rows);
    kernel.setArg(3, src_offset_cols);
    kernel.setArg(4, dst_offset_rows);
    kernel.setArg(5, dst_offset_cols);

    kernel.setArg(6, size_rows);
    kernel.setArg(7, size_cols);

    kernel.setArg(8, src.rows());
    kernel.setArg(9, src.cols());
    kernel.setArg(10, dst.rows());
    kernel.setArg(11, dst.cols());

    cmdQueue.enqueueNDRangeKernel(kernel, cl::NullRange,
                                  cl::NDRange(size_rows, size_cols),
                                  cl::NullRange, NULL, NULL);
  } catch (const cl::Error& e) {
    check_opencl_error("copy_submatrix", e);
  }
}
}  // namespace math
}  // namespace stan

#endif
#endif
