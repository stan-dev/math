#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <complex>

// Diagonal shifts keep lambda_i(A) + mu_j(B) well away from zero, so the
// problems are well conditioned. Random nonsymmetric matrices still have
// complex eigenvalue pairs, so the complex Schur path is exercised.

TEST(MathFunctions, solve_sylvester_real_known_solution) {
  using stan::math::solve_sylvester;
  const int n = 6;
  const int m = 5;
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(n, n)
                      + 2.0 * n * Eigen::MatrixXd::Identity(n, n);
  Eigen::MatrixXd B = Eigen::MatrixXd::Random(m, m)
                      + 2.0 * m * Eigen::MatrixXd::Identity(m, m);
  Eigen::MatrixXd X_true = Eigen::MatrixXd::Random(n, m);
  Eigen::MatrixXd C = A * X_true + X_true * B;

  Eigen::MatrixXcd X = solve_sylvester(A, B, C);

  EXPECT_EQ(X.rows(), n);
  EXPECT_EQ(X.cols(), m);
  // real inputs give a real solution up to roundoff
  EXPECT_LT(X.imag().norm(), 1e-12);
  EXPECT_LT((X.real() - X_true).norm(), 1e-10);
  EXPECT_LT((A * X + X * B - C).norm(), 1e-10);
}

TEST(MathFunctions, solve_sylvester_complex_known_solution) {
  using stan::math::solve_sylvester;
  const int n = 6;
  const int m = 5;
  Eigen::MatrixXcd A = Eigen::MatrixXcd::Random(n, n)
                       + 2.0 * n * Eigen::MatrixXcd::Identity(n, n);
  Eigen::MatrixXcd B = Eigen::MatrixXcd::Random(m, m)
                       + 2.0 * m * Eigen::MatrixXcd::Identity(m, m);
  Eigen::MatrixXcd X_true = Eigen::MatrixXcd::Random(n, m);
  Eigen::MatrixXcd C = A * X_true + X_true * B;

  Eigen::MatrixXcd X = solve_sylvester(A, B, C);

  EXPECT_LT((X - X_true).norm(), 1e-10);
  EXPECT_LT((A * X + X * B - C).norm(), 1e-10);
}

TEST(MathFunctions, solve_sylvester_mixed_real_complex) {
  using stan::math::solve_sylvester;
  const int n = 4;
  const int m = 3;
  Eigen::MatrixXd A_r = Eigen::MatrixXd::Random(n, n)
                        + 2.0 * n * Eigen::MatrixXd::Identity(n, n);
  Eigen::MatrixXd B_r = Eigen::MatrixXd::Random(m, m)
                        + 2.0 * m * Eigen::MatrixXd::Identity(m, m);
  Eigen::MatrixXcd A_c = Eigen::MatrixXcd::Random(n, n)
                         + 2.0 * n * Eigen::MatrixXcd::Identity(n, n);

  // real A and B, complex C: known solution
  Eigen::MatrixXcd X_true = Eigen::MatrixXcd::Random(n, m);
  Eigen::MatrixXcd C_c = A_r * X_true + X_true * B_r;
  Eigen::MatrixXcd X1 = solve_sylvester(A_r, B_r, C_c);
  EXPECT_LT((X1 - X_true).norm(), 1e-10);

  // complex A, real B and C: residual
  Eigen::MatrixXd C_r = Eigen::MatrixXd::Random(n, m);
  Eigen::MatrixXcd X2 = solve_sylvester(A_c, B_r, C_r);
  EXPECT_LT((A_c * X2 + X2 * B_r - C_r).norm(), 1e-10);
}

TEST(MathFunctions, solve_sylvester_diagonal) {
  // for diagonal A and B: X(i, j) = C(i, j) / (A(i, i) + B(j, j))
  Eigen::VectorXd a(3);
  a << 1, 2, 3;
  Eigen::VectorXd b(2);
  b << 4, 5;
  Eigen::MatrixXd A = a.asDiagonal();
  Eigen::MatrixXd B = b.asDiagonal();
  Eigen::MatrixXd C(3, 2);
  C << 1, 2, 3, 4, 5, 6;

  Eigen::MatrixXcd X = stan::math::solve_sylvester(A, B, C);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 2; ++j) {
      EXPECT_NEAR(X(i, j).real(), C(i, j) / (a(i) + b(j)), 1e-12);
      EXPECT_NEAR(X(i, j).imag(), 0.0, 1e-12);
    }
  }
}

TEST(MathFunctions, solve_sylvester_expression_inputs) {
  // Eigen expressions go through plain_type_t; result must match plain inputs
  using stan::math::solve_sylvester;
  Eigen::MatrixXd A
      = Eigen::MatrixXd::Random(4, 4) + 8.0 * Eigen::MatrixXd::Identity(4, 4);
  Eigen::MatrixXd B
      = Eigen::MatrixXd::Random(3, 3) + 6.0 * Eigen::MatrixXd::Identity(3, 3);
  Eigen::MatrixXd C = Eigen::MatrixXd::Random(4, 3);
  Eigen::MatrixXd At = A.transpose();

  Eigen::MatrixXcd X_expr = solve_sylvester(A.transpose(), B, 2.0 * C);
  Eigen::MatrixXcd X_plain = solve_sylvester(At, B, Eigen::MatrixXd(2.0 * C));

  EXPECT_LT((X_expr - X_plain).norm(), 1e-14);
}

TEST(MathFunctions, solve_sylvester_errors) {
  using stan::math::solve_sylvester;
  Eigen::MatrixXd A = Eigen::MatrixXd::Identity(3, 3);
  Eigen::MatrixXd B = Eigen::MatrixXd::Identity(2, 2);
  Eigen::MatrixXd C = Eigen::MatrixXd::Ones(3, 2);

  // valid input does not throw
  EXPECT_NO_THROW(solve_sylvester(A, B, C));

  // A not square
  Eigen::MatrixXd A_rect = Eigen::MatrixXd::Ones(3, 2);
  EXPECT_THROW(solve_sylvester(A_rect, B, C), std::invalid_argument);

  // B not square
  Eigen::MatrixXd B_rect = Eigen::MatrixXd::Ones(2, 3);
  EXPECT_THROW(solve_sylvester(A, B_rect, C), std::invalid_argument);

  // rows of C do not match A
  Eigen::MatrixXd C_rows = Eigen::MatrixXd::Ones(4, 2);
  EXPECT_THROW(solve_sylvester(A, B, C_rows), std::invalid_argument);

  // cols of C do not match B
  Eigen::MatrixXd C_cols = Eigen::MatrixXd::Ones(3, 3);
  EXPECT_THROW(solve_sylvester(A, B, C_cols), std::invalid_argument);
}

TEST(MathFunctions, solve_sylvester_zero_size) {
  using stan::math::solve_sylvester;
  Eigen::MatrixXd I2 = Eigen::MatrixXd::Identity(2, 2);
  Eigen::MatrixXd E(0, 0);

  Eigen::MatrixXcd X1 = solve_sylvester(E, I2, Eigen::MatrixXd(0, 2));
  EXPECT_EQ(X1.rows(), 0);
  EXPECT_EQ(X1.cols(), 2);

  Eigen::MatrixXcd X2 = solve_sylvester(I2, E, Eigen::MatrixXd(2, 0));
  EXPECT_EQ(X2.rows(), 2);
  EXPECT_EQ(X2.cols(), 0);
}

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
