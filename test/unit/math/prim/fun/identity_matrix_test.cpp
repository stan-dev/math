#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, identity_matrix) {
  using Eigen::MatrixXd;
  using stan::math::identity_matrix;

  MatrixXd u0 = identity_matrix(0);
  EXPECT_EQ(0, u0.rows());
  EXPECT_EQ(0, u0.cols());

  for (int K = 1; K < 5; K++) {
    MatrixXd m = Eigen::MatrixXd::Identity(K, K);
    EXPECT_MATRIX_FLOAT_EQ(m, identity_matrix(K));
  }
}

TEST(MathFunctions, identity_matrix_throw) {
  EXPECT_THROW(stan::math::identity_matrix(-1), std::domain_error);
}
