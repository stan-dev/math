#ifndef STAN_MATH_GPU_CHOLESKY_DECOMPOSE_HPP
#define STAN_MATH_GPU_CHOLESKY_DECOMPOSE_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/matrix_gpu.hpp>
#include <stan/math/gpu/kernels/cholesky_decompose.hpp>
#include <stan/math/gpu/multiply.hpp>
#include <stan/math/gpu/multiply_transpose.hpp>
#include <stan/math/gpu/lower_tri_inverse.hpp>
#include <stan/math/gpu/transpose.hpp>
#include <stan/math/gpu/subtract.hpp>
#include <stan/math/gpu/err/check_diagonal_zeros.hpp>
#include <stan/math/gpu/err/check_nan.hpp>
#include <CL/cl.hpp>
#include <algorithm>
#include <cmath>

namespace stan {
namespace math {
inline matrix_gpu cholesky_decompose(matrix_gpu& A, const int min_block);
namespace internal {
inline matrix_gpu cholesky_decompose_recursion(matrix_gpu& A,
                                               const int min_block) {
  matrix_gpu L(A.rows(), A.cols());
  if (A.rows() <= min_block || A.rows() < 100) {
    try {
      opencl_kernels::cholesky_decompose(cl::NDRange(A.rows()),
                                         cl::NDRange(A.rows()), A.buffer(),
                                         L.buffer(), A.rows());
    } catch (const cl::Error& e) {
      check_opencl_error("cholesky_decompose", e);
    }
  } else {
    L = stan::math::cholesky_decompose(A, min_block);
  }
  return L;
}
}  // namespace internal
/**
 * Return the lower-triangular Cholesky factor (i.e., matrix
 * square root) of the specified square, symmetric matrix.
 * The return value \f$L\f$ will be a lower-traingular matrix such that the
 * original matrix \f$A\f$ is given by
 * <p>\f$A = L \times L^T\f$.
 * The Cholesky decomposition is computed on the GPU. This algorithm is
 * recursive. The matrix is subset into a matrix of size
 *  <code>A.rows() / 2</code>, and if the <code>block</code> size is less than
 * 100 then the cholesky decomposition on the GPU is computed
 * using that submatrix. If <code>block</code> is greater than
 * 100 then <code>cholesky_decompose</code> is run again
 * with <code>block</code> equal to <code>A.rows() / 2</code>. Once the
 * Cholesky Decomposition is computed, the full matrix cholesky is created
 * by propogating the cholesky forward as given in the reference report below.
 *
 * For a full guide to how this works
 * see the Cholesy decompostion chapter in the  reference report
 * <a href="https://goo.gl/6kWkJ5"> here</a>.
 * @param A Symmetric matrix on the GPU.
 * @return Square root of matrix on the GPU.
 * @throw std::domain_error if m is not
 *  positive definite (if m has more than 0 elements)
 */
inline matrix_gpu cholesky_decompose(matrix_gpu& A, const int min_block) {
  auto offset = 0;
  auto block = floor(A.rows() / 2);
  // NOTE: The code in this section follows the naming conventions
  // in the report linked in the docs.
  matrix_gpu A_11(block, block);
  // Repeats the blocked cholesky decomposition until the size of the remaining
  // submatrix is smaller or equal to the block size
  if ((offset + block) < (A.rows())) {
    auto block_subset = A.rows() - offset - block;
    matrix_gpu A_21(block_subset, block);
    matrix_gpu A_22(block_subset, block_subset);
    // Copies a block of the input A into A_11
    A_11.sub_block(A, offset, offset, 0, 0, block, block);
    // The following function either calls the
    // blocked cholesky recursively for the submatrix A_11
    // or calls the kernel  directly if the size of the block is small enough
    matrix_gpu L_11
        = stan::math::internal::cholesky_decompose_recursion(A_11, min_block);
    // Copies L_11 back to the input matrix
    A.sub_block(L_11, 0, 0, offset, offset, block, block);
    // Copies a block of the input A into A_21
    auto block_offset = offset + block;
    A_21.sub_block(A, block_offset, offset, 0, 0, block_subset, block);
    // computes A_21*((L_11^-1)^T)
    // and copies the resulting submatrix to the input matrix
    matrix_gpu L_21 = A_21 * transpose(lower_triangular_inverse(L_11));
    A.sub_block(L_21, 0, 0, block_offset, offset, block_subset, block);
    A_22.sub_block(A, block_offset, block_offset, 0, 0, block_subset,
                   block_subset);
    // computes A_22 - L_21*(L_21^T)
    matrix_gpu L_22 = A_22 - multiply_transpose(L_21);
    A.sub_block(L_22, 0, 0, block_offset, block_offset, block_subset,
                block_subset);
    offset += block;
  }
  // Computes the Cholesky factor for the remaining submatrix
  const auto remaining_rows = A.rows() - offset;
  if (remaining_rows > 0) {
    matrix_gpu A_11(remaining_rows, remaining_rows);
    A_11.sub_block(A, offset, offset, 0, 0, remaining_rows, remaining_rows);
    // calculate the cholesky factor for the remaining part of the matrix
    matrix_gpu L_11
        = stan::math::internal::cholesky_decompose_recursion(A_11, min_block);
    A.sub_block(L_11, 0, 0, offset, offset, remaining_rows, remaining_rows);
  }
  check_nan("cholesky_decompose_gpu", "Matrix m", A);
  check_diagonal_zeros("cholesky_decompose_gpu", "Matrix m", A);
  A.zeros<stan::math::TriangularViewGPU::Upper>();
  return A;
}
}  // namespace math
}  // namespace stan

#endif
#endif
