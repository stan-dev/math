#include <stan/math/prim/mat/fun/mrrr.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, tridiag_eigensolver_trivial) {
  Eigen::VectorXd diag(3), subdiag(2), eigenvals;
  diag << 1.5, 1.2, -2;
  subdiag << 0, 0;

  Eigen::MatrixXd eigenvecs(3, 3);
  stan::math::internal::tridiagonal_eigensolver(diag, subdiag, eigenvals,
                                                eigenvecs);
  EXPECT_TRUE(eigenvals.isApprox(diag));
  EXPECT_TRUE(eigenvecs.isApprox(Eigen::MatrixXd::Identity(3, 3)));
}

TEST(MathMatrix, tridiag_eigensolver_small) {
  int size = 7;
  srand(0);  // ensure test repeatability
  Eigen::VectorXd diag = Eigen::VectorXd::Random(size);
  Eigen::VectorXd subdiag = Eigen::VectorXd::Random(size - 1);
  subdiag[2] = 0;

  Eigen::VectorXd eigenvals;
  Eigen::MatrixXd eigenvecs;
  stan::math::internal::tridiagonal_eigensolver(diag, subdiag, eigenvals,
                                                eigenvecs);

  Eigen::MatrixXd t = Eigen::MatrixXd::Constant(size, size, 0);
  t.diagonal() = diag;
  t.diagonal(1) = subdiag;
  t.diagonal(-1) = subdiag;

  EXPECT_NEAR(diag.sum(), eigenvals.sum(), 1e-11);
  EXPECT_TRUE((eigenvecs * eigenvecs.transpose())
                  .isApprox(Eigen::MatrixXd::Identity(size, size)));
  EXPECT_TRUE((t * eigenvecs).isApprox(eigenvecs * eigenvals.asDiagonal()));
}

TEST(MathMatrix, tridiag_eigensolver_large) {
  int size = 345;
  srand(0);  // ensure test repeatability
  Eigen::VectorXd diag = Eigen::VectorXd::Random(size);
  Eigen::VectorXd subdiag = Eigen::VectorXd::Random(size - 1);
  subdiag[12] = 0;
  subdiag[120] = 0;
  subdiag[121] = 0;

  Eigen::VectorXd eigenvals;
  Eigen::MatrixXd eigenvecs;
  stan::math::internal::tridiagonal_eigensolver(diag, subdiag, eigenvals,
                                                eigenvecs);

  Eigen::MatrixXd t = Eigen::MatrixXd::Constant(size, size, 0);
  t.diagonal() = diag;
  t.diagonal(1) = subdiag;
  t.diagonal(-1) = subdiag;

  EXPECT_NEAR(diag.sum(), eigenvals.sum(), 1e-9);
  EXPECT_TRUE((eigenvecs * eigenvecs.transpose())
                  .isApprox(Eigen::MatrixXd::Identity(size, size), 1e-9));
  EXPECT_TRUE(
      (t * eigenvecs).isApprox(eigenvecs * eigenvals.asDiagonal(), 1e-9));
}
