#ifndef STAN_MATH_OPENCL_SUB_BLOCK_HPP
#define STAN_MATH_OPENCL_SUB_BLOCK_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernels/sub_block.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <cl.hpp>
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
template <typename T>
inline void matrix_cl<T, require_arithmetic_t<T>>::sub_block(
    const matrix_cl<T, require_arithmetic_t<T>>& A, size_t A_i, size_t A_j,
    size_t this_i, size_t this_j, size_t nrows, size_t ncols) try {
  if (nrows == 0 || ncols == 0) {
    return;
  }
  if ((A_i + nrows) > A.rows() || (A_j + ncols) > A.cols()
      || (this_i + nrows) > this->rows() || (this_j + ncols) > this->cols()) {
    domain_error("sub_block", "submatrix in *this", " is out of bounds", "");
  }
  cl::CommandQueue cmdQueue = opencl_context.queue();
  if (A.view() == matrix_cl_view::Entire) {
    std::array<size_t, 3> src_offset({A_i * sizeof(double), A_j, 0});
    std::array<size_t, 3> dst_offset({this_i * sizeof(double), this_j, 0});
    std::array<size_t, 3> size({nrows * sizeof(double), ncols, 1});
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
                              this->rows(), this->cols(), A.view());
  }
  // calculation of extreme sub- and super- diagonal written
  const int diag_in_copy = A_i - A_j;
  const int copy_low = contains_nonzero(A.view(), matrix_cl_view::Lower)
                           ? 1 - nrows
                           : diag_in_copy;
  const int copy_high = contains_nonzero(A.view(), matrix_cl_view::Upper)
                            ? ncols - 1
                            : diag_in_copy;
  const int start = this_j - this_i;

  if (start + copy_low < 0) {
    this->view_ = either(this->view_, matrix_cl_view::Lower);
  } else if (this_i <= 1 && this_j == 0 && nrows + this_i >= rows_
             && ncols >= std::min(rows_, cols_) - 1
             && !contains_nonzero(A.view_, matrix_cl_view::Lower)) {
    this->view_ = both(this->view_, matrix_cl_view::Upper);
  }

  if (start + copy_high > 0) {
    this->view_ = either(this->view_, matrix_cl_view::Upper);
  } else if (this_i == 0 && this_j <= 1 && ncols + this_j >= cols_
             && nrows >= std::min(rows_, cols_) - 1
             && !contains_nonzero(A.view_, matrix_cl_view::Upper)) {
    this->view_ = both(this->view_, matrix_cl_view::Lower);
  }
} catch (const cl::Error& e) {
  check_opencl_error("copy_submatrix", e);
}

}  // namespace math
}  // namespace stan

#endif
#endif
