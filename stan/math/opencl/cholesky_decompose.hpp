#ifndef STAN_MATH_OPENCL_CHOLESKY_DECOMPOSE_HPP
#define STAN_MATH_OPENCL_CHOLESKY_DECOMPOSE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/err.hpp>
#include <stan/math/opencl/multiply_transpose.hpp>
#include <stan/math/opencl/prim/multiply.hpp>
#include <stan/math/opencl/tri_inverse.hpp>
#include <stan/math/opencl/kernels/cholesky_decompose.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/meta.hpp>
#include <CL/opencl.hpp>
#include <algorithm>
#include <cmath>

namespace stan {
namespace math {
namespace opencl {
/** \ingroup opencl
 * Performs an in-place computation of the lower-triangular Cholesky factor
 * (i.e., matrix square root) of the specified square, symmetric matrix. The
 * return value \f$L\f$ will be a lower-triangular matrix such that the original
 * matrix \f$A\f$ is given by <p>\f$A = L \times L^T\f$. The Cholesky
 * decomposition is computed using an OpenCL kernel. This algorithm is
 * recursive. The matrix is subset into a matrix of size <code>A.rows() /
 * 4</code>, and if the submatrix size is less than 50 or <code>min_block</code>
 * then the Cholesky decomposition on the OpenCL device is computed using that
 * submatrix. If the submatrix is greater than 50 or <code>min_block</code> then
 * <code>cholesky_decompose</code> is run again on a submatrix with size equal
 * to <code>submat.rows() / 4</code>. Once the Cholesky decomposition is
 * computed, the full matrix Cholesky is created by propagating the Cholesky
 * forward as given in the reference report below.
 *
 * For a full guide to how this works
 * see the Cholesky decomposition chapter in the reference report
 * <a href="https://goo.gl/6kWkJ5"> here</a>.
 * @throw std::domain_error if m is not
 *  positive definite (if m has more than 0 elements)
 */
template <typename T, typename = require_floating_point_t<T>>
inline void cholesky_decompose(matrix_cl<T>& A) {
  if (A.rows() == 0) {
    return;
  }
  // Repeats the blocked cholesky decomposition until the size of the remaining
  // submatrix is smaller or equal to the minimum blocks size
  // or a heuristic of 100.
  // The Cholesky (OpenCL) algorithm only uses one local block so we need the
  // matrix To be less than the max thread block size.
  if (A.rows() <= opencl_context.tuning_opts().cholesky_min_L11_size) {
    try {
      opencl_kernels::cholesky_decompose(cl::NDRange(A.rows()),
                                         cl::NDRange(A.rows()), A, A.rows());
    } catch (const cl::Error& e) {
      check_opencl_error("cholesky_decompose", e);
    }
    A.view(matrix_cl_view::Lower);
    return;
  }
  // NOTE: The code in this section follows the naming conventions
  // in the report linked in the docs.
  const int block
      = std::floor(A.rows() / opencl_context.tuning_opts().cholesky_partition);
  // Subset the top left block of the input A into A_11
  matrix_cl<T> A_11 = block_zero_based(A, 0, 0, block, block);
  // The following function either calls the
  // blocked cholesky recursively for the submatrix A_11
  // or calls the kernel  directly if the size of the block is small enough
  opencl::cholesky_decompose(A_11);
  // Copies L_11 back to the input matrix
  block_zero_based(A, 0, 0, block, block)
      = block_zero_based(A_11, 0, 0, block, block);

  const int block_subset = A.rows() - block;
  matrix_cl<T> A_21 = block_zero_based(A, block, 0, block_subset, block);
  // computes A_21*((L_11^-1)^T)
  // and copies the resulting submatrix to the lower left hand corner of A
  matrix_cl<T> L_21 = A_21 * transpose(tri_inverse(A_11));
  block_zero_based(A, block, 0, block_subset, block)
      = block_zero_based(L_21, 0, 0, block_subset, block);
  matrix_cl<T> A_22
      = block_zero_based(A, block, block, block_subset, block_subset);
  // computes A_22 - L_21*(L_21^T)
  matrix_cl<T> L_22 = A_22 - multiply_transpose(L_21);
  // copy L_22 into A's lower left hand corner
  opencl::cholesky_decompose(L_22);
  block_zero_based(A, block, block, block_subset, block_subset)
      = block_zero_based(L_22, 0, 0, block_subset, block_subset);
  A.view(matrix_cl_view::Lower);
}
}  // namespace opencl
}  // namespace math
}  // namespace stan

#endif
#endif
