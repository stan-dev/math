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
inline void block_householder_tridiag_cl(const matrix_cl<double>& A,
                                         matrix_cl<double>& packed,
                                         const int r = 60) {
  packed = A;
  for (size_t k = 0; k < A.rows() - 2; k += r) {
    const int actual_r = std::min({r, static_cast<int>(A.rows() - k - 2)});
    matrix_cl<double> V_cl = constant(0.0, A.rows() - k - 1, actual_r + 1);

    matrix_cl<double> Uu(actual_r, 1), Vu(actual_r, 1), q_cl(1, 1);
    for (size_t j = 0; j < actual_r; j++) {
      try {
        int hh_local
            = opencl_kernels::tridiagonalization_householder.get_option(
                "LOCAL_SIZE_");
        opencl_kernels::tridiagonalization_householder(
            cl::NDRange(hh_local), cl::NDRange(hh_local), packed, V_cl, q_cl,
            packed.rows(), V_cl.rows(), j, k);
        if (j != 0) {
          int v_step_1_local
              = opencl_kernels::tridiagonalization_v_step_1.get_option(
                  "LOCAL_SIZE_");
          opencl_kernels::tridiagonalization_v_step_1(
              cl::NDRange(v_step_1_local * j), cl::NDRange(v_step_1_local),
              packed, V_cl, Uu, Vu, packed.rows(), V_cl.rows(), k);
        }
        int v_step_2_local
            = opencl_kernels::tridiagonalization_v_step_2.get_option(
                "LOCAL_SIZE_");
        opencl_kernels::tridiagonalization_v_step_2(
            cl::NDRange((A.rows() - k - j - 1 + v_step_2_local - 1)
                        / v_step_2_local * v_step_2_local),
            cl::NDRange(v_step_2_local), packed, V_cl, Uu, Vu, packed.rows(),
            V_cl.rows(), k, j);
        int v_step_3_local
            = opencl_kernels::tridiagonalization_v_step_3.get_option(
                "LOCAL_SIZE_");
        opencl_kernels::tridiagonalization_v_step_3(
            cl::NDRange(v_step_3_local), cl::NDRange(v_step_3_local), packed,
            V_cl, q_cl, packed.rows(), V_cl.rows(), k, j);
      } catch (cl::Error& e) {
        check_opencl_error("block_householder_tridiag_cl", e);
      }
    }
    matrix_cl<double> U_cl = block_zero_based(
        packed, k + actual_r, k, A.rows() - k - actual_r, actual_r);
    matrix_cl<double> V_block_cl = block_zero_based(
        V_cl, actual_r - 1, 0, V_cl.rows() - actual_r + 1, actual_r);
    matrix_cl<double> partial_update_cl = U_cl * transpose(V_block_cl);

    auto block
        = block_zero_based(packed, k + actual_r, k + actual_r,
                           partial_update_cl.rows(), partial_update_cl.cols());
    block = block - partial_update_cl - transpose(partial_update_cl);
  }
  block_zero_based(packed, packed.rows() - 2, packed.cols() - 1, 1, 1)
      = block_zero_based(packed, packed.rows() - 1, packed.cols() - 2, 1, 1);
}

/**
 * Calculates Q*A in-place. To construct Q pass an appropriate identity matrix
 * as input A.
 * @param packed Packed result of tridiagonalization that contains householder
 * vectors that define Q in columns bellow the diagonal. Usually result of a
 * call to `block_householder_tridiag_cl`.
 * @param[in,out] A On input a matrix to multiply with Q. On output the product
 * Q*A.
 * @param r Block size. Affects only performance of the algorithm. Optimal value
 * depends on the size of A and cache of the processor. For larger matrices or
 * larger cache sizes a larger value is optimal.
 */
inline void block_apply_packed_Q_cl(const matrix_cl<double>& packed_cl,
                                    matrix_cl<double>& A, const int r = 200) {
  Eigen::MatrixXd packed = from_matrix_cl(packed_cl);
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
    matrix_cl<double> packed_block_transpose_triang_cl = transpose(
        block_zero_based(packed_cl, k + 1, k, packed.rows() - k - 1, actual_r));
    packed_block_transpose_triang_cl.view(matrix_cl_view::Upper);
    matrix_cl<double> W_cl(W);
    auto A_bottom_cl
        = block_zero_based(A, k + 1, 0, A.rows() - k - 1, A.cols());
    matrix_cl<double> A_bottom_cl_eval = A_bottom_cl;
    matrix_cl<double> tmp1
        = packed_block_transpose_triang_cl * A_bottom_cl_eval;
    matrix_cl<double> tmp2 = W_cl * tmp1;
    A_bottom_cl -= tmp2;
  }
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
#endif
