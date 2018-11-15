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
namespace stan {
namespace math {
/**
 * Return the lower-triangular Cholesky factor (i.e., matrix
 * square root) of the specified square, symmetric matrix.
 * The return value \f$L\f$ will be a lower-traingular matrix such that the
 * original matrix \f$A\f$ is given by
 * <p>\f$A = L \times L^T\f$.
 * The Cholesky decomposition is computed on the GPU. This algorithm is
 * recursive, where The parameters <code>block</code>, <code>divider</code>, and
 *  <code>min_block</code> act as tuning parameters for the recursive step of
 *  the GPU based Cholesky decompostion. The matrix is subset by the
 *  <code>block</code> size, and if the <code>block</code> size is less than
 * <code>min_block</code> then the cholesky decomposition on the GPU is computed
 * using that submatrix. If <code>block</code> is greater than
 * <code>block_size</code> then <code>cholesky_decompose</code> is run again
 * with <code>block</code> equal to <code>block/divider</code>. Once the
 * Cholesky Decomposition is computed, the full matrix cholesky is created
 * by propogating the cholesky forward as given in the reference report below.
 *
 * For a full guide to how this works
 * see the Cholesy decompostion chapter in the  reference report
 * <a href="https://goo.gl/6kWkJ5"> here</a>.
 * @param A Symmetric matrix on the GPU.
 * @param block Size of the block used to compute the cholesky decomposition.
 * @param divider Proportion to divide the submatrix by at each recursive step.
 * @param min_block The amount that block is checked against to decide
 *  whether to continue the recursion or to perform the cholesky.
 * @return Square root of matrix on the GPU.
 * @throw std::domain_error if m is not
 *  positive definite (if m has more than 0 elements)
 */
inline matrix_gpu cholesky_decompose(matrix_gpu& A, const int block = 100,
                                     const int divider = 2,
                                     const int min_block = 100) {
  auto offset = 0;
  // NOTE: The code in this section follows the naming conventions
  // in the report linked in the docs.
  matrix_gpu A_11(block, block);
  matrix_gpu L_11(block, block);
  // Repeats the blocked cholesky decomposition until the size of the remaining
  // submatrix is smaller or equal to the block size
  while ((offset + block) < (A.rows())) {
    auto block_subset = A.rows() - offset - block;
    matrix_gpu A_21(block_subset, block);
    matrix_gpu A_22(block_subset, block_subset);
    // Copies a block of the input A into A_11
    A_11.sub_block(A, offset, offset, 0, 0, block, block);
    // Calls the blocked cholesky for the submatrix A_11
    // or calls the kernel  directly if the size of the block is small enough
    if (block <= min_block || divider <= 1) {
      try {
        opencl_kernels::cholesky_decompose(
            cl::NDRange(A_11.rows()), cl::NDRange(A_11.rows()), A_11.buffer(),
            L_11.buffer(), A_11.rows());
      } catch (const cl::Error& e) {
        check_opencl_error("cholesky_decompose", e);
      }
    } else {
      L_11 = cholesky_decompose(A_11, block / divider, divider, min_block);
    }
    // Copies L_11 back to the input matrix
    A.sub_block(L_11, 0, 0, offset, offset, block, block);
    // Copies a block of the input A into A_21
    auto block_offset = offset + block;
    A_21.sub_block(A, block_offset, offset, 0, 0, block_subset, block);
    // computes A_21*((L_11^-1)^T)
    // and copies the resulting submatrix to the input matrix
    auto A_11_inverse_T = transpose(lower_triangular_inverse(L_11));
    // TODO(Steve): Replace with mult operator when that PR goes through
    matrix_gpu L_21 = multiply(A_21, A_11_inverse_T);
    A.sub_block(L_21, 0, 0, block_offset, offset, block_subset, block);
    A_22.sub_block(A, block_offset, block_offset, 0, 0, block_subset,
                   block_subset);
    // computes A_22 - L_21*(L_21^T)
    matrix_gpu temp = multiply_transpose(L_21);
    // TODO(Steve): Replace with subtraction operator when that PR goes through
    matrix_gpu L_22 = subtract(A_22, temp);
    A.sub_block(L_22, 0, 0, block_offset, block_offset, block_subset,
                block_subset);
    offset += block;
  }
  // Computes the Cholesky factor for the remaining part of the matrix
  const auto remaining_rows = A.rows() - offset;
  if (remaining_rows > 0) {
    matrix_gpu A_11(remaining_rows, remaining_rows);
    matrix_gpu L_11(remaining_rows, remaining_rows);
    A_11.sub_block(A, offset, offset, 0, 0, remaining_rows, remaining_rows);
    // Calls the blocked cholesky for the submatrix A_11
    // or calls the kernel  directly if the size of the block is small enough
    if (block <= min_block || divider <= 1) {
      try {
        opencl_kernels::cholesky_decompose(
            cl::NDRange(A_11.rows()), cl::NDRange(A_11.rows()), A_11.buffer(),
            L_11.buffer(), A_11.rows());
      } catch (const cl::Error& e) {
        check_opencl_error("cholesky_decompose", e);
      }
    } else {
      L_11 = cholesky_decompose(A_11, block / divider);
    }
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
