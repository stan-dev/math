#ifndef STAN_MATH_OPENCL_TRI_INVERSE_HPP
#define STAN_MATH_OPENCL_TRI_INVERSE_HPP

#ifdef STAN_OPENCL

#include <stan/math/opencl/plain_type.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernels/diag_inv.hpp>
#include <stan/math/opencl/kernels/inv_lower_tri_multiply.hpp>
#include <stan/math/opencl/kernels/neg_rect_lower_tri_multiply.hpp>
#include <stan/math/opencl/err.hpp>
#include <stan/math/opencl/kernels/batch_identity.hpp>
#include <stan/math/opencl/zeros_strict_tri.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/meta.hpp>
#include <cmath>
#include <string>
#include <vector>

namespace stan {
namespace math {
/** \ingroup opencl
 * Computes the inverse of a triangular matrix
 *
 * For a full guide to how this works and fits into Cholesky decompositions,
 * see the reference report
 * <a href="https://github.com/SteveBronder/stancon2018/blob/master/report.pdf">
 * here</a> and kernel doc
 * <a href="https://github.com/stan-dev/math/wiki/GPU-Kernels">here</a>.
 *
 * @param A matrix on the OpenCL device
 * @return the inverse of A
 *
 * @throw <code>std::invalid_argument</code> if the matrix
 *    is not square
 */
template <matrix_cl_view matrix_view = matrix_cl_view::Entire, typename T,
          require_matrix_cl_st<std::is_floating_point, T>* = nullptr>
inline plain_type_t<T> tri_inverse(const T& A) {
  check_square("tri_inverse (OpenCL)", "A", A);
  // if the triangular view is not specified use the triangularity of
  // the input matrix
  matrix_cl_view tri_view = matrix_view;
  if (matrix_view == matrix_cl_view::Entire) {
    if (A.view() != matrix_cl_view::Diagonal) {
      check_triangular("tri_inverse (OpenCL)", "A", A);
    }
    tri_view = A.view();
  }
  if (tri_view == matrix_cl_view::Diagonal) {
    plain_type_t<T> inv_mat(A.rows(), A.cols());
    diagonal(inv_mat) = elt_divide(1.0, diagonal(A));
    return inv_mat;
  }

  int thread_block_2D_dim = 32;
  int max_1D_thread_block_size = opencl_context.max_thread_block_size();
  // we split the input matrix to 32 blocks
  int thread_block_size_1D
      = (((A.rows() / 32) + thread_block_2D_dim - 1) / thread_block_2D_dim)
        * thread_block_2D_dim;
  if (max_1D_thread_block_size < thread_block_size_1D) {
    thread_block_size_1D = max_1D_thread_block_size;
  }
  int max_2D_thread_block_dim = std::sqrt(max_1D_thread_block_size);
  if (max_2D_thread_block_dim < thread_block_2D_dim) {
    thread_block_2D_dim = max_2D_thread_block_dim;
  }
  // for small size split in max 2 parts
  if (thread_block_size_1D < 64) {
    thread_block_size_1D = 32;
  }
  if (A.rows() < thread_block_size_1D) {
    thread_block_size_1D = A.rows();
  }

  // pad the input matrix
  int A_rows_padded
      = ((A.rows() + thread_block_size_1D - 1) / thread_block_size_1D)
        * thread_block_size_1D;

  plain_type_t<T> temp(A_rows_padded, A_rows_padded);
  plain_type_t<T> inv_padded = constant(0.0, A_rows_padded, A_rows_padded);
  plain_type_t<T> inv_mat(A);
  plain_type_t<T> zero_mat
      = constant(0.0, A_rows_padded - A.rows(), A_rows_padded);
  if (tri_view == matrix_cl_view::Upper) {
    inv_mat = transpose(inv_mat).eval();
  }
  int work_per_thread
      = opencl_kernels::inv_lower_tri_multiply.get_option("WORK_PER_THREAD");
  // the number of blocks in the first step
  // each block is inverted with using the regular forward substitution
  int parts = inv_padded.rows() / thread_block_size_1D;
  block_zero_based(inv_padded, 0, 0, inv_mat.rows(), inv_mat.rows()) = inv_mat;
  try {
    // create a batch of identity matrices to be used in the first step
    opencl_kernels::batch_identity(
        cl::NDRange(parts, thread_block_size_1D, thread_block_size_1D), temp,
        thread_block_size_1D, temp.size());
    // spawn parts thread blocks, each responsible for one block
    opencl_kernels::diag_inv(cl::NDRange(parts * thread_block_size_1D),
                             cl::NDRange(thread_block_size_1D), inv_padded,
                             temp, inv_padded.rows());
  } catch (cl::Error& e) {
    check_opencl_error("inverse step1", e);
  }
  // set the padded part of the matrix and the upper triangular to zeros
  block_zero_based(inv_padded, inv_mat.rows(), 0, zero_mat.rows(),
                   zero_mat.cols())
      = zero_mat;
  inv_padded.template zeros_strict_tri<stan::math::matrix_cl_view::Upper>();
  if (parts == 1) {
    inv_mat
        = block_zero_based(inv_padded, 0, 0, inv_mat.rows(), inv_mat.rows());
    if (tri_view == matrix_cl_view::Upper) {
      inv_mat = transpose(inv_mat).eval();
    }
    return inv_mat;
  }
  using std::ceil;
  parts = ceil(parts / 2.0);

  auto result_matrix_dim = thread_block_size_1D;
  auto thread_block_work2d_dim = thread_block_2D_dim / work_per_thread;
  auto ndrange_2d
      = cl::NDRange(thread_block_2D_dim, thread_block_work2d_dim, 1);
  while (parts > 0) {
    int result_matrix_dim_x = result_matrix_dim;
    // when calculating the last submatrix
    // we can reduce the size to the actual size (not the next power of 2)
    if (parts == 1 && (inv_padded.rows() - result_matrix_dim * 2) < 0) {
      result_matrix_dim_x = inv_padded.rows() - result_matrix_dim;
    }
    auto result_work_dim = result_matrix_dim / work_per_thread;
    auto result_ndrange
        = cl::NDRange(result_matrix_dim_x, result_work_dim, parts);
    opencl_kernels::inv_lower_tri_multiply(result_ndrange, ndrange_2d,
                                           inv_padded, temp, inv_padded.rows(),
                                           result_matrix_dim);
    opencl_kernels::neg_rect_lower_tri_multiply(
        result_ndrange, ndrange_2d, inv_padded, temp, inv_padded.rows(),
        result_matrix_dim);
    // if this is the last submatrix, end
    if (parts == 1) {
      parts = 0;
    } else {
      parts = ceil(parts / 2.0);
    }
    result_matrix_dim *= 2;
    // set the padded part and upper diagonal to zeros
    block_zero_based(inv_padded, inv_mat.rows(), 0, zero_mat.rows(),
                     zero_mat.cols())
        = zero_mat;
    inv_padded.template zeros_strict_tri<stan::math::matrix_cl_view::Upper>();
  }
  // un-pad and return
  inv_mat = block_zero_based(inv_padded, 0, 0, inv_mat.rows(), inv_mat.rows());
  if (tri_view == matrix_cl_view::Upper) {
    inv_mat = transpose(inv_mat).eval();
  }
  inv_mat.view(tri_view);
  return inv_mat;
}
}  // namespace math
}  // namespace stan

#endif
#endif
