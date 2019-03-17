#ifndef STAN_MATH_OPENCL_SUB_BLOCK_HPP
#define STAN_MATH_OPENCL_SUB_BLOCK_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/constants.hpp>
#include <stan/math/opencl/kernels/sub_block.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <CL/cl.hpp>

namespace stan {
namespace math {

/**
 * Write the context of A into
 * <code>this</code> starting at the top left of <code>this</code>
 * @param A input matrix
 * @param A_i the offset row in A
 * @param A_j the offset column in A
 * @param this_i the offset row for the matrix to be subset into
 * @param this_j the offset col for the matrix to be subset into
 * @param nrows the number of rows in the submatrix
 * @param ncols the number of columns in the submatrix
 */
template <TriangularViewCL triangular_view = TriangularViewCL::Entire>
inline void matrix_cl::sub_block(const matrix_cl& A, int A_i, int A_j, int this_i, int this_j,
               int nrows, int ncols) {
  if (nrows == 0 || ncols == 0) {
    return;
  }
  if ((A_i + nrows) > A.rows() || (A_j + ncols) > A.cols()
      || (this_i + nrows) > this->rows() || (this_j + ncols) > this->cols()) {
    domain_error("sub_block", "submatrix in *this", " is out of bounds", "");
  }
  try {
    cl::Event sub_block_event = opencl_kernels::sub_block(cl::NDRange(nrows, ncols), A,
                              *this, A_i, A_j, this_i, this_j, nrows,
                              ncols, A.rows(), A.cols(), this->rows(),
                              this->cols(), triangular_view);
    this->events(sub_block_event);
  } catch (const cl::Error& e) {
    check_opencl_error("copy_submatrix", e);
  }
}

}  // namespace math
}  // namespace stan

#endif
#endif
