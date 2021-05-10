#ifndef STAN_MATH_OPENCL_TRIDIAGONALIZATION_HPP
#define STAN_MATH_OPENCL_TRIDIAGONALIZATION_HPP

#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/prim/multiply.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/copy.hpp>

#include <stan/math/opencl/kernels/tridiagonalization.hpp>

namespace stan {
namespace math {
namespace internal {

/**
 * Tridiagonalize a symmetric matrix using block Housholder algorithm. A = Q * T
 * * Q^T, where T is tridiagonal and Q is orthonormal.
 * @param A Input matrix
 * @param[out] packed Packed form of the tridiagonal matrix. Elements of the
 * resulting symmetric tridiagonal matrix T are in the diagonal and first
 * superdiagonal. Columns bellow diagonal contain householder vectors that can
 * be used to construct orthogonal matrix Q.
 * @param r Block size. Affects only performance of the algorithm. Optimal value
 * depends on the size of A and cache of the processor. For larger matrices or
 * larger cache sizes a larger value is optimal.
 */
void block_householder_tridiag_cl(const Eigen::MatrixXd& A,
                                  Eigen::MatrixXd& packed, const int r = 60) {
  matrix_cl<double> packed_gpu(A);
  for (size_t k = 0; k < A.rows() - 2; k += r) {
    const int actual_r = std::min({r, static_cast<int>(A.rows() - k - 2)});
    matrix_cl<double> V_gpu(A.rows() - k - 1, actual_r + 1);

    matrix_cl<double> Uu(actual_r, 1), Vu(actual_r, 1), q_gpu(1, 1);
    for (size_t j = 0; j < actual_r; j++) {
      try {
        int hh_local
            = opencl_kernels::tridiagonalization_householder.get_option(
                "LOCAL_SIZE_");
        opencl_kernels::tridiagonalization_householder(
            cl::NDRange(hh_local), cl::NDRange(hh_local), packed_gpu, V_gpu,
            q_gpu, packed_gpu.rows(), V_gpu.rows(), j, k);
        if (j != 0) {
          int v1_local
              = opencl_kernels::tridiagonalization_v1.get_option("LOCAL_SIZE_");
          opencl_kernels::tridiagonalization_v1(
              cl::NDRange(v1_local * j), cl::NDRange(v1_local), packed_gpu,
              V_gpu, Uu, Vu, packed_gpu.rows(), V_gpu.rows(), k);
        }
        int v2_local
            = opencl_kernels::tridiagonalization_v2.get_option("LOCAL_SIZE_");
        opencl_kernels::tridiagonalization_v2(
            cl::NDRange((A.rows() - k - j - 1 + v2_local - 1) / v2_local
                        * v2_local),
            cl::NDRange(v2_local), packed_gpu, V_gpu, Uu, Vu, packed_gpu.rows(),
            V_gpu.rows(), k, j);
        int v3_local
            = opencl_kernels::tridiagonalization_v3.get_option("LOCAL_SIZE_");
        opencl_kernels::tridiagonalization_v3(
            cl::NDRange(v3_local), cl::NDRange(v3_local), packed_gpu, V_gpu,
            q_gpu, packed_gpu.rows(), V_gpu.rows(), k, j);
      } catch (cl::Error& e) {
        check_opencl_error("block_householder_tridiag_cl", e);
      }
    }
    matrix_cl<double> U_gpu = block_zero_based(
        packed_gpu, k + actual_r, k, A.rows() - k - actual_r, actual_r);
    matrix_cl<double> V_block_gpu = block_zero_based(
        V_gpu, actual_r - 1, 0, V_gpu.rows() - actual_r + 1, actual_r);
    matrix_cl<double> partial_update_gpu = U_gpu * transpose(V_block_gpu);

    auto block = block_zero_based(packed_gpu, k + actual_r, k + actual_r,
                                  partial_update_gpu.rows(),
                                  partial_update_gpu.cols());
    block = block - partial_update_gpu - transpose(partial_update_gpu);
  }
  packed = from_matrix_cl(packed_gpu);
  packed(packed.rows() - 2, packed.cols() - 1)
      = packed(packed.rows() - 1, packed.cols() - 2);
}

/**
 * Calculates Q*A in place. To construct Q pass an appropriate identity matrix
 * as input A.
 * @param packed Packed result of tridiagonalization that contains householder
 * vectors that define Q in columns bellow the diagonal. Usually result of a
 * call to `block_householder_tridiag_cl`.
 * @param[in,out] A On input a matrix to multiply with Q. On output the product
 * Q*A.
 * @param r Block size. Affects only performance of the algorithm. Optimal value
 * depends on the size of A and cache of the processor. For larger matrices or
 * larger cache sizes larger value is optimal.
 */
void block_apply_packed_Q_cl(const Eigen::MatrixXd& packed, Eigen::MatrixXd& A,
                             const int r = 200) {
  matrix_cl<double> A_gpu(A);
  Eigen::MatrixXd scratch_space(A.rows(), r);
  for (int k = (packed.rows() - 3) / r * r; k >= 0; k -= r) {
    const int actual_r = std::min({r, static_cast<int>(packed.rows() - k - 2)});
    Eigen::MatrixXd W(packed.rows() - k - 1, actual_r);
    W.col(0) = packed.col(k).tail(W.rows());
    for (size_t j = 1; j < actual_r; j++) {
      scratch_space.col(0).head(j).noalias()
          = packed.block(k + j + 1, k, packed.rows() - k - j - 1, j).transpose()
            * packed.col(j + k).tail(packed.rows() - k - j - 1);
      W.col(j).noalias() = -W.leftCols(j) * scratch_space.col(0).head(j);
      W.col(j).tail(W.rows() - j)
          += packed.col(j + k).tail(packed.rows() - k - j - 1);
    }
    Eigen::MatrixXd packed_block_transpose_triang
        = packed.block(k + 1, k, packed.rows() - k - 1, actual_r)
              .transpose()
              .triangularView<Eigen::Upper>();
    matrix_cl<double> packed_block_transpose_triang_gpu(
        packed_block_transpose_triang, matrix_cl_view::Upper);
    matrix_cl<double> W_gpu(W);
    auto A_bottom_gpu
        = block_zero_based(A_gpu, k + 1, 0, A.rows() - k - 1, A.cols());
    matrix_cl<double> A_bottom_gpu_eval = A_bottom_gpu;
    matrix_cl<double> tmp1 = packed_block_transpose_triang_gpu * A_bottom_gpu_eval;
    matrix_cl<double> tmp2 = W_gpu * tmp1;
    A_bottom_gpu -= tmp2;
  }
  A = from_matrix_cl(A_gpu);
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
#endif
