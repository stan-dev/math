#ifndef STAN_MATH_PRIM_MAT_FUN_TRIDIAGONALIZATION_HPP
#define STAN_MATH_PRIM_MAT_FUN_TRIDIAGONALIZATION_HPP

#include <Eigen/Dense>

namespace stan {
namespace math {

/**
 * Tridiagonalize a symmetric matrix using block Housholder algorithm. A = Q * T * Q^T, where T is tridiagonal and Q is orthonormal.
 * @param A Input matrix
 * @param[out] packed Packed form of the tridiagonal matrix. Elements of the resulting symmetric tridiagonal matrix T are in the diagonal and first superdiagonal.
 * Columns bellow diagonal contain householder vectors that can be used to construct orthogonal matrix Q.
 * @param r Block size. Affects only performance of the algorithm. Optimal value depends on the size of A and cache of the processor. For larger matrices or larger cache sizes a larger value is optimal.
 */
void block_householder_tridiag(const Eigen::MatrixXd& A, Eigen::MatrixXd& packed, int r = 60) {
  packed = A;
  for (size_t k = 0; k < packed.rows() - 2; k += r) {
    int actual_r = std::min({r, static_cast<int>(packed.rows() - k - 2)});
    Eigen::MatrixXd V(packed.rows() - k - 1, actual_r);

    for (size_t j = 0; j < actual_r; j++) {
      auto householder = packed.col(k + j).tail(packed.rows() - k - j - 1);
      if (j != 0) {
        auto householder_whole = packed.col(k + j).tail(packed.rows() - k - j);
        householder_whole -= packed.block(j + k, k, householder_whole.size(), j) * V.block(j - 1, 0, 1, j).transpose() +
                             V.block(j - 1, 0, householder_whole.size(), j) * packed.block(j + k, k, 1, j).transpose();
      }
      double q = householder.squaredNorm();
      double alpha = -copysign(sqrt(q), packed(k + j, k + j));
      q -= householder[0] * householder[0];
      householder[0] -= alpha;
      q += householder[0] * householder[0];
      q = sqrt(q);
      householder *= SQRT_2 / q;

      auto& u = householder;
      Eigen::VectorXd v(householder.size() + 1);
      v.tail(householder.size()) = packed.bottomRightCorner(packed.rows() - k - j - 1, packed.cols() - k - j - 1).selfadjointView<Eigen::Lower>() * u
                                   - packed.block(k + j + 1, k, u.size(), j) * (V.bottomLeftCorner(u.size(), j).transpose() * u)
                                   - V.bottomLeftCorner(u.size(), j) * (packed.block(k + j + 1, k, u.size(), j).transpose() * u);
      v[0] = q / SQRT_2;
      double cnst = v.tail(householder.size()).transpose() * u;
      v.tail(householder.size()) -= 0.5 * cnst * u;

      packed(k + j, k + j + 1) = packed(k + j + 1, k + j) * q / SQRT_2 + alpha - v[0] * u[0];
      V.col(j).tail(V.rows() - j) = v.tail(V.rows() - j);
    }
    Eigen::MatrixXd partial_update = packed.block(k + actual_r, k,packed.rows() - k - actual_r, actual_r) * V.bottomRows(V.rows() - actual_r + 1).transpose();
    packed.block(k + actual_r, k + actual_r, packed.rows() - k - actual_r, packed.cols() - k - actual_r).triangularView<Eigen::Lower>() -= partial_update + partial_update.transpose();
  }
  packed(packed.rows() - 2, packed.cols() - 1) = packed(packed.rows() - 1, packed.cols() - 2);
}

/**
 * Calculates Q*A in place. To construct Q pass identity matrix as input A.
 * @param packed Packed result of tridiagonalization that contains householder vectors that define Q in columns bellow the diagonal. Usually result of a call to `block_householder_tridiag3`.
 * @param[in,out] A On input a matrix to multiply with Q. On output the product Q*A.
 * @param r Block size. Affects only performance of the algorithm. Optimal value depends on the size of A and cache of the processor. For larger matrices or larger cache sizes a larger value is optimal.
 */
void block_apply_packed_Q(const Eigen::MatrixXd& packed, Eigen::MatrixXd& A, int r = 100) {
  //if input A==Identity, constructs Q
  Eigen::MatrixXd scratchSpace(A.rows(), r);
  for (int k = (packed.rows() - 3) / r * r; k >= 0; k -= r) {
    int actual_r = std::min({r, static_cast<int>(packed.rows() - k - 2)});
    Eigen::MatrixXd W(packed.rows() - k - 1, actual_r);
    W.col(0) = packed.col(k).tail(W.rows());
    for (size_t j = 1; j < actual_r; j++) {
      scratchSpace.col(0).head(j).noalias() = packed.block(k + j + 1, k, packed.rows() - k - j - 1, j).transpose() * packed.col(j + k).tail(packed.rows() - k - j - 1);
      W.col(j).noalias() = -W.leftCols(j) * scratchSpace.col(0).head(j);
      W.col(j).tail(W.rows() - j) += packed.col(j + k).tail(packed.rows() - k - j - 1);
    }
    scratchSpace.transpose().bottomRows(actual_r).noalias() = packed.block(k + 1, k, packed.rows() - k - 1, actual_r).transpose().triangularView<Eigen::Upper>() * A.bottomRows(A.rows() - k - 1);
    A.bottomRows(A.cols() - k - 1).noalias() -= W * scratchSpace.transpose().bottomRows(actual_r);
  }
}

}
}

#endif
