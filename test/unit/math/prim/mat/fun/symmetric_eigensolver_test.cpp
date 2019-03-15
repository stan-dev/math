#include <stan/math/prim/mat/fun/symmetric_eigensolver.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, symmetric_eigensolver_trivial) {
  Eigen::MatrixXd input(3, 3);
  input << 3, 0, 0, 0, 2, 0, 0, 0, 1;

  Eigen::VectorXd eigenvals;
  Eigen::MatrixXd eigenvecs;
  stan::math::symmetric_eigensolver(input, eigenvals, eigenvecs);

  EXPECT_TRUE(eigenvals.isApprox(input.diagonal()));
  EXPECT_TRUE(eigenvecs.isApprox(Eigen::MatrixXd::Identity(3, 3)));
}

TEST(MathMatrix, symmetric_eigensolver_small) {
  int size = 7;
  srand(0);  // ensure test repeatability
  Eigen::MatrixXd input = Eigen::MatrixXd::Random(size, size);
  input += input.transpose().eval();

  Eigen::VectorXd eigenvals;
  Eigen::MatrixXd eigenvecs;
  stan::math::symmetric_eigensolver(input, eigenvals, eigenvecs);

  EXPECT_NEAR(input.diagonal().sum(), eigenvals.sum(), 1e-11);
  EXPECT_TRUE((eigenvecs * eigenvecs.transpose())
                  .isApprox(Eigen::MatrixXd::Identity(size, size)));
  EXPECT_TRUE((input * eigenvecs).isApprox(eigenvecs * eigenvals.asDiagonal()));
}

TEST(MathMatrix, symmetric_eigensolver_large) {
  int size = 345;
  srand(0);  // ensure test repeatability
  Eigen::MatrixXd input = Eigen::MatrixXd::Random(size, size);
  input += input.transpose().eval();

  Eigen::VectorXd eigenvals;
  Eigen::MatrixXd eigenvecs;
  stan::math::symmetric_eigensolver(input, eigenvals, eigenvecs);

  EXPECT_NEAR(input.diagonal().sum(), eigenvals.sum(), 1e-9);
  EXPECT_TRUE((eigenvecs * eigenvecs.transpose())
                  .isApprox(Eigen::MatrixXd::Identity(size, size), 1e-9));
  EXPECT_TRUE((input * eigenvecs).isApprox(eigenvecs * eigenvals.asDiagonal()));
}
