#ifndef STAN_MATH_OPENCL_QR_DECOMPOSITION_HPP
#define STAN_MATH_OPENCL_QR_DECOMPOSITION_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/zeros_strict_tri.hpp>
#include <stan/math/opencl/prim/identity_matrix.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/multiply.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/opencl/copy.hpp>
#include <CL/opencl.hpp>
#include <algorithm>
#include <vector>
#include <cmath>

namespace stan {
namespace math {

/**
 * Calculates QR decomposition of A using the block Householder algorithm.
 *
 * The implementation is based on the article:
 * \link
 * https://cpb-us-w2.wpmucdn.com/sites.gatech.edu/dist/5/462/files/2016/08/Kerr_Campbell_Richards_QRD_on_GPUs.pdf
 * @tparam need_Q whether Q needs to be calculated (default true)
 * @param A matrix to factorize
 * @param Q out the orthonormal matrix
 * @param R out the lower triangular matrix
 * @param r Block size. Optimal value depends on the hardware.
 */
template <bool need_Q = true>
void qr_decomposition_cl(const matrix_cl<double>& A, matrix_cl<double>& Q,
                         matrix_cl<double>& R, int r = 100) {
  using std::copysign;
  using std::sqrt;
  int rows = A.rows();
  int cols = A.cols();
  if (need_Q) {
    Q = identity_matrix<matrix_cl<double>>(rows);
  }
  R = A;
  if (rows <= 1) {
    return;
  }

  matrix_cl<double> V_alloc(rows, r);
  matrix_cl<double> W_alloc(rows, r);
  Eigen::VectorXd temporary(rows > cols ? rows : cols);

  for (int k = 0; k < std::min(rows - 1, cols); k += r) {
    int actual_r = std::min({r, rows - k - 1, cols - k});
    matrix_cl<double> V(V_alloc.buffer(), rows - k, actual_r);
    matrix_cl<double>& Y_cl = V;
    matrix_cl<double> W_cl(W_alloc.buffer(), rows - k, actual_r);

    Eigen::MatrixXd R_cpu_block
        = from_matrix_cl(block_zero_based(R, k, k, rows - k, actual_r));

    Eigen::MatrixXd V_cpu(rows - k, actual_r);
    V_cpu.triangularView<Eigen::StrictlyUpper>()
        = Eigen::MatrixXd::Constant(V.rows(), V.cols(), 0);

    for (int j = 0; j < actual_r; j++) {
      Eigen::VectorXd householder = R_cpu_block.block(j, j, rows - k - j, 1);
      householder.coeffRef(0)
          -= copysign(householder.norm(), householder.coeff(0));
      if (householder.rows() != 1) {
        householder *= SQRT_TWO / householder.norm();
      }

      V_cpu.col(j).tail(V.rows() - j) = householder;

      Eigen::VectorXd::AlignedMapType temp_map(temporary.data(), actual_r - j);
      temp_map.noalias()
          = R_cpu_block.block(j, j, rows - k - j, actual_r - j).transpose()
            * householder;
      R_cpu_block.block(j, j + 1, rows - k - j, actual_r - j - 1).noalias()
          -= householder * temp_map.tail(actual_r - j - 1).transpose();
      R_cpu_block.coeffRef(j, j) -= householder.coeff(0) * temp_map.coeff(0);
    }

    matrix_cl<double> R_back_block = to_matrix_cl(R_cpu_block);
    R_back_block.view(matrix_cl_view::Upper);
    block_zero_based(R, k, k, rows - k, actual_r) = R_back_block;

    Eigen::MatrixXd& Y_cpu = V_cpu;
    Eigen::MatrixXd W_cpu = V_cpu;
    for (int j = 1; j < actual_r; j++) {
      Eigen::VectorXd::AlignedMapType temp_map(temporary.data(), j);
      temp_map.noalias() = Y_cpu.leftCols(j).transpose() * V_cpu.col(j);
      W_cpu.col(j).noalias() -= W_cpu.leftCols(j) * temp_map;
    }
    W_cl = to_matrix_cl(W_cpu);
    V = to_matrix_cl(V_cpu);

    auto R_block
        = block_zero_based(R, k, k + actual_r, rows - k, cols - k - actual_r);
    R_block -= Y_cl * (transpose(W_cl) * R_block);

    if (need_Q) {
      auto Q_block = block_zero_based(Q, 0, k, Q.rows(), Q.cols() - k);
      Q_block -= (Q_block * W_cl) * transpose(Y_cl);
    }
  }

  R.view(matrix_cl_view::Upper);
}

}  // namespace math
}  // namespace stan
#endif
#endif  // STAN_MATH_OPENCL_QR_DECOMPOSITION_HPP
