#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

namespace {

// sparse block-diagonal matrix with the block pattern the Laplace solver
// reserves for W_r (see CholeskyWSolverBlock)
Eigen::SparseMatrix<double> block_pattern(int n_blocks, int block_size) {
  const int n = n_blocks * block_size;
  Eigen::SparseMatrix<double> m(n, n);
  m.reserve(Eigen::VectorXi::Constant(n, block_size));
  for (int b = 0; b < n_blocks; ++b) {
    for (int k = 0; k < block_size; ++k) {
      for (int j = 0; j < block_size; ++j) {
        m.insert(b * block_size + j, b * block_size + k) = 1.0;
      }
    }
  }
  m.makeCompressed();
  return m;
}

Eigen::SparseMatrix<double> block_diag(
    const std::vector<Eigen::MatrixXd>& blocks) {
  const int block_size = blocks[0].rows();
  Eigen::SparseMatrix<double> w = block_pattern(blocks.size(), block_size);
  for (std::size_t b = 0; b < blocks.size(); ++b) {
    for (int k = 0; k < block_size; ++k) {
      for (int j = 0; j < block_size; ++j) {
        w.coeffRef(b * block_size + j, b * block_size + k) = blocks[b](j, k);
      }
    }
  }
  return w;
}

void expect_principal_sqrt(const Eigen::SparseMatrix<double>& w,
                           int block_size) {
  Eigen::SparseMatrix<double> w_root = block_pattern(w.rows() / block_size,
                                                     block_size);
  EXPECT_NO_THROW(
      stan::math::internal::block_matrix_sqrt(w_root, w, block_size));
  const Eigen::MatrixXd root = w_root;
  EXPECT_TRUE(root.isApprox(root.transpose(), 1e-12));
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(root);
  EXPECT_GE(eig.eigenvalues().minCoeff(), -1e-12);
  EXPECT_TRUE((root * root).isApprox(Eigen::MatrixXd(w), 1e-10));
}

}  // namespace

TEST(LaplaceBlockMatrixSqrt, PositiveDefiniteBlocks) {
  Eigen::MatrixXd m(3, 3);
  m << 1.0, 0.3, -0.2, 0.5, 2.0, 0.1, -0.4, 0.2, 1.5;
  Eigen::MatrixXd a = m * m.transpose() + Eigen::MatrixXd::Identity(3, 3);
  Eigen::MatrixXd b = 2.0 * Eigen::MatrixXd::Identity(3, 3);
  expect_principal_sqrt(block_diag({a, b}), 3);
}

// A negative Hessian with more latent variables than observations is only
// positive semi-definite: its zero eigenvalues come out of floating point
// as tiny values of either sign and must not be rejected.
TEST(LaplaceBlockMatrixSqrt, RankDeficientBlockIsAccepted) {
  Eigen::MatrixXd z(3, 6);
  z << 1.0, 0.5, -0.3, 0.8, 0.1, -0.6, 0.2, -1.1, 0.4, 0.3, 0.9, 0.7, -0.5, 0.6,
      1.2, -0.2, 0.4, 0.3;
  Eigen::MatrixXd w = z.transpose() * z;  // rank 3 of 6
  expect_principal_sqrt(block_diag({w}), 6);

  Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(2, 2);  // rank 1 of 2
  Eigen::MatrixXd spd = Eigen::MatrixXd::Identity(2, 2);
  expect_principal_sqrt(block_diag({spd, ones}), 2);
}

TEST(LaplaceBlockMatrixSqrt, IndefiniteBlockThrows) {
  Eigen::MatrixXd indefinite(2, 2);
  indefinite << 1.0, 0.0, 0.0, -1.0;
  Eigen::SparseMatrix<double> w = block_diag({indefinite});
  Eigen::SparseMatrix<double> w_root = block_pattern(1, 2);
  EXPECT_THROW(stan::math::internal::block_matrix_sqrt(w_root, w, 2),
               std::domain_error);
}

TEST(LaplaceBlockMatrixSqrt, NonFiniteBlockThrows) {
  Eigen::MatrixXd nan_block = Eigen::MatrixXd::Identity(2, 2);
  nan_block(0, 1) = std::numeric_limits<double>::quiet_NaN();
  Eigen::SparseMatrix<double> w = block_diag({nan_block});
  Eigen::SparseMatrix<double> w_root = block_pattern(1, 2);
  EXPECT_THROW(stan::math::internal::block_matrix_sqrt(w_root, w, 2),
               std::domain_error);
}
