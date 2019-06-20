#ifndef STAN_MATH_OPENCL_SUB_BLOCK_HPP
#define STAN_MATH_OPENCL_SUB_BLOCK_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/triangular.hpp>
#include <stan/math/opencl/kernels/sub_block.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <CL/cl.hpp>
#include <vector>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Write the contents of A into
 * <code>this</code> starting at the top left of <code>this</code>
 * @param A input matrix
 * @param A_i the offset row in A
 * @param A_j the offset column in A
 * @param this_i the offset row for the matrix to be subset into
 * @param this_j the offset col for the matrix to be subset into
 * @param nrows the number of rows in the submatrix
 * @param ncols the number of columns in the submatrix
 */
inline void matrix_cl::sub_block(const matrix_cl& A, size_t A_i, size_t A_j,
                                 size_t this_i, size_t this_j, size_t nrows,
                                 size_t ncols) try {
  if (nrows == 0 || ncols == 0) {
    return;
  }
  if ((A_i + nrows) > A.rows() || (A_j + ncols) > A.cols()
      || (this_i + nrows) > this->rows() || (this_j + ncols) > this->cols()) {
    domain_error("sub_block", "submatrix in *this", " is out of bounds", "");
  }
  cl::CommandQueue cmdQueue = opencl_context.queue();
  if (A.triangular_view() == TriangularViewCL::Entire) {
    cl::size_t<3> src_offset
        = opencl::to_size_t<3>({A_i * sizeof(double), A_j, 0});
    cl::size_t<3> dst_offset
        = opencl::to_size_t<3>({this_i * sizeof(double), this_j, 0});
    cl::size_t<3> size
        = opencl::to_size_t<3>({nrows * sizeof(double), ncols, 1});
    std::vector<cl::Event> kernel_events
        = vec_concat(A.write_events(), this->read_write_events());
    cl::Event copy_event;
    cmdQueue.enqueueCopyBufferRect(A.buffer(), this->buffer(), src_offset,
                                   dst_offset, size, A.rows() * sizeof(double),
                                   A.rows() * A.cols() * sizeof(double),
                                   sizeof(double) * this->rows(),
                                   this->rows() * this->cols() * sizeof(double),
                                   &kernel_events, &copy_event);
    A.add_read_event(copy_event);
    this->add_write_event(copy_event);
  } else {
    opencl_kernels::sub_block(cl::NDRange(nrows, ncols), A, *this, A_i, A_j,
                              this_i, this_j, nrows, ncols, A.rows(), A.cols(),
                              this->rows(), this->cols(), A.triangular_view());
  }
  // calculation of extreme sub- and super- diagonal written
  int diag_in_copy = A_i - A_j;
  int copy_low = static_cast<bool>(A.triangular_view_ & TriangularViewCL::Lower)
                     ? 1 - nrows
                     : diag_in_copy;
  int copy_high
      = static_cast<bool>(A.triangular_view_ & TriangularViewCL::Upper)
            ? ncols - 1
            : diag_in_copy;
  int start = this_j - this_i;

  if (start + copy_low < 0) {
    triangular_view_ |= TriangularViewCL::Lower;
  } else if (this_i <= 1 && this_j == 0 && nrows + this_i >= rows_
             && ncols >= std::min(rows_, cols_) - 1
             && !static_cast<bool>(A.triangular_view_
                                   & TriangularViewCL::Lower)) {
    triangular_view_ &= TriangularViewCL::Upper;
  }

  if (start + copy_high > 0) {
    triangular_view_ |= TriangularViewCL::Upper;
  } else if (this_i == 0 && this_j <= 1 && ncols + this_j >= cols_
             && nrows >= std::min(rows_, cols_) - 1
             && !static_cast<bool>(A.triangular_view_
                                   & TriangularViewCL::Upper)) {
    triangular_view_ &= TriangularViewCL::Lower;
  }
} catch (const cl::Error& e) {
  check_opencl_error("copy_submatrix", e);
}

}  // namespace math
}  // namespace stan

#endif
#endif
