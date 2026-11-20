#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <complex>

// Stable A (by Gershgorin all eigenvalues have Re <= -(n + 1) < 0) keeps
// lambda_i + conj(lambda_j) away from zero.

TEST(MathFunctions, solve_lyapunov_real_known_solution) {
  using stan::math::solve_lyapunov;
  const int n = 6;
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(n, n)
                      - 2.0 * n * Eigen::MatrixXd::Identity(n, n);
  Eigen::MatrixXd X_true = Eigen::MatrixXd::Random(n, n);
  Eigen::MatrixXd C = A * X_true + X_true * A.transpose();

  Eigen::MatrixXcd X = solve_lyapunov(A, C);

  EXPECT_EQ(X.rows(), n);
  EXPECT_EQ(X.cols(), n);
  // real inputs give a real solution up to roundoff
  EXPECT_LT(X.imag().norm(), 1e-12);
  EXPECT_LT((X.real() - X_true).norm(), 1e-10);
  EXPECT_LT((A * X + X * A.transpose() - C).norm(), 1e-10);
}

TEST(MathFunctions, solve_lyapunov_complex_known_solution) {
  using stan::math::solve_lyapunov;
  const int n = 6;
  Eigen::MatrixXcd A = Eigen::MatrixXcd::Random(n, n)
                       - 2.0 * n * Eigen::MatrixXcd::Identity(n, n);
  Eigen::MatrixXcd X_true = Eigen::MatrixXcd::Random(n, n);
  Eigen::MatrixXcd C = A * X_true + X_true * A.adjoint();

  Eigen::MatrixXcd X = solve_lyapunov(A, C);

  EXPECT_LT((X - X_true).norm(), 1e-10);
  EXPECT_LT((A * X + X * A.adjoint() - C).norm(), 1e-10);
}

TEST(MathFunctions, solve_lyapunov_matches_sylvester) {
  using stan::math::solve_lyapunov;
  using stan::math::solve_sylvester;
  const int n = 5;
  Eigen::MatrixXcd A = Eigen::MatrixXcd::Random(n, n)
                       - 2.0 * n * Eigen::MatrixXcd::Identity(n, n);
  Eigen::MatrixXcd C = Eigen::MatrixXcd::Random(n, n);
  Eigen::MatrixXcd A_adj = A.adjoint();

  Eigen::MatrixXcd X_lyap = solve_lyapunov(A, C);
  Eigen::MatrixXcd X_sylv = solve_sylvester(A, A_adj, C);

  EXPECT_LT((X_lyap - X_sylv).norm(), 1e-12);

  // real A, complex C
  Eigen::MatrixXd A_r = A.real();
  Eigen::MatrixXd A_rt = A_r.transpose();
  EXPECT_LT((solve_lyapunov(A_r, C) - solve_sylvester(A_r, A_rt, C)).norm(),
            1e-12);
}

TEST(MathFunctions, solve_lyapunov_hermitian) {
  // stable A and Hermitian negative definite C = -Q give a Hermitian
  // positive definite solution
  using stan::math::solve_lyapunov;
  const int n = 5;
  Eigen::MatrixXcd A = Eigen::MatrixXcd::Random(n, n)
                       - 2.0 * n * Eigen::MatrixXcd::Identity(n, n);
  Eigen::MatrixXcd L = Eigen::MatrixXcd::Random(n, n);
  Eigen::MatrixXcd Q = L * L.adjoint() + Eigen::MatrixXcd::Identity(n, n);

  Eigen::MatrixXcd X = solve_lyapunov(A, -Q);

  EXPECT_LT((A * X + X * A.adjoint() + Q).norm(), 1e-10);
  EXPECT_LT((X - X.adjoint()).norm(), 1e-12);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eig(0.5 * (X + X.adjoint()));
  EXPECT_GT(eig.eigenvalues().minCoeff(), 0.0);
}

TEST(MathFunctions, solve_lyapunov_errors) {
  using stan::math::solve_lyapunov;
  Eigen::MatrixXd A = Eigen::MatrixXd::Identity(3, 3);
  Eigen::MatrixXd C = Eigen::MatrixXd::Ones(3, 3);

  // valid input does not throw
  EXPECT_NO_THROW(solve_lyapunov(A, C));

  // A not square
  Eigen::MatrixXd A_rect = Eigen::MatrixXd::Ones(3, 2);
  EXPECT_THROW(solve_lyapunov(A_rect, C), std::invalid_argument);

  // rows of C do not match A
  Eigen::MatrixXd C_rows = Eigen::MatrixXd::Ones(4, 3);
  EXPECT_THROW(solve_lyapunov(A, C_rows), std::invalid_argument);

  // cols of C do not match A
  Eigen::MatrixXd C_cols = Eigen::MatrixXd::Ones(3, 4);
  EXPECT_THROW(solve_lyapunov(A, C_cols), std::invalid_argument);
}

TEST(MathFunctions, solve_lyapunov_zero_size) {
  Eigen::MatrixXd E(0, 0);
  Eigen::MatrixXcd X = stan::math::solve_lyapunov(E, E);
  EXPECT_EQ(X.rows(), 0);
  EXPECT_EQ(X.cols(), 0);
}
