#ifndef STAN_MATH_GPU_LOWER_TRI_INVERSE_HPP
#define STAN_MATH_GPU_LOWER_TRI_INVERSE_HPP

#ifdef STAN_OPENCL
#include <stan/math/gpu/matrix_gpu.hpp>
#include <stan/math/gpu/kernels/lower_tri_inverse_step1.hpp>
#include <stan/math/gpu/kernels/lower_tri_inverse_step2.hpp>
#include <stan/math/gpu/kernels/lower_tri_inverse_step3.hpp>

#include <stan/math/gpu/identity.hpp>
#include <stan/math/gpu/err/check_square.hpp>
#include <string>
#include <vector>

namespace stan {
namespace math {
/**
 * Computes the inverse of the lower triangular matrix
 * that resides in the GPU global memory
 *
 * @param A matrix on the GPU
 *
 * @return the inverse of A
 *
 * @throw <code>std::invalid_argument</code> if the matrix
 *    is not square
 */
inline matrix_gpu lower_triangular_inverse(const matrix_gpu& A) {
  check_square("lower_triangular_inverse (GPU)", "A", A);

  int thread_block_size_1D = 32;
  int thread_block_2D_dim = 32;
  int A_rows_padded
      = ((A.rows() + thread_block_2D_dim - 1) / thread_block_2D_dim)
        * thread_block_2D_dim;

  matrix_gpu temp(A_rows_padded, A_rows_padded);
  matrix_gpu inv_padded(A_rows_padded, A_rows_padded);
  matrix_gpu inv(A);
  inv_padded.zeros<stan::math::TriangularViewGPU::Entire>();

  int max_1D_thread_block_size = opencl_context.max_thread_block_size();
  int max_2D_thread_block_dim = sqrt(max_1D_thread_block_size);
  if (max_1D_thread_block_size < thread_block_size_1D) {
    thread_block_size_1D = max_1D_thread_block_size;
  }
  if (inv_padded.rows() < thread_block_size_1D) {
    thread_block_size_1D = inv_padded.rows();
  }
  if (max_2D_thread_block_dim < thread_block_2D_dim) {
    thread_block_2D_dim = max_2D_thread_block_dim;
  }

  int wpt = opencl_kernels::lower_tri_inverse_step2.make_functor.get_opts().at(
      "WORK_PER_THREAD");
  int parts
      = (inv_padded.rows() + thread_block_size_1D - 1) / thread_block_size_1D;
  try {
    opencl_kernels::lower_tri_inverse_step1(
        cl::NDRange(parts * thread_block_size_1D),
        cl::NDRange(thread_block_size_1D), inv.buffer(), temp.buffer(),
        inv.rows());
  } catch (cl::Error& e) {
    check_opencl_error("inverse step1", e);
  }
  parts /= 2;
  if (parts == 0) {
    return inv;
  }
  inv_padded.sub_block(inv, 0, 0, 0, 0, inv.rows(), inv.rows());

  int result_matrix_dim = thread_block_size_1D;
  while (parts > 1) {
    opencl_kernels::lower_tri_inverse_step2(
        cl::NDRange(result_matrix_dim, result_matrix_dim / wpt, parts),
        cl::NDRange(thread_block_2D_dim, thread_block_2D_dim / wpt, 1),
        inv_padded.buffer(), temp.buffer(), inv_padded.rows(),
        result_matrix_dim, result_matrix_dim);
    opencl_kernels::lower_tri_inverse_step3(
        cl::NDRange(result_matrix_dim, result_matrix_dim / wpt, parts),
        cl::NDRange(thread_block_2D_dim, thread_block_2D_dim / wpt, 1),
        inv_padded.buffer(), temp.buffer(), inv_padded.rows(),
        result_matrix_dim, result_matrix_dim);
    parts /= 2;
    result_matrix_dim *= 2;
  }
  opencl_kernels::lower_tri_inverse_step2(
      cl::NDRange(inv.rows() - result_matrix_dim, result_matrix_dim / wpt,
                  parts),
      cl::NDRange(thread_block_2D_dim, thread_block_2D_dim / wpt, 1),
      inv_padded.buffer(), temp.buffer(), inv_padded.rows(),
      inv.rows() - result_matrix_dim, result_matrix_dim);
  opencl_kernels::lower_tri_inverse_step3(
      cl::NDRange(inv.rows() - result_matrix_dim, result_matrix_dim / wpt,
                  parts),
      cl::NDRange(thread_block_2D_dim, thread_block_2D_dim / wpt, 1),
      inv_padded.buffer(), temp.buffer(), inv_padded.rows(), result_matrix_dim,
      result_matrix_dim);
  inv.sub_block(inv_padded, 0, 0, 0, 0, A.rows(), A.rows());
  return inv;
}
}  // namespace math
}  // namespace stan

#endif
#endif
