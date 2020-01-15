#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, identity_matrix) {
  using Eigen::MatrixXd;
  using stan::math::identity_matrix;

  MatrixXd u0 = identity_matrix(0);
  EXPECT_EQ(0, u0.rows());
  EXPECT_EQ(0, u0.cols());

  MatrixXd u1 = identity_matrix(1);
  MatrixXd m1(1, 1);
  m1 << 1;

  EXPECT_EQ(1, u1.rows());
  EXPECT_EQ(1, u1.cols());
  EXPECT_EQ(u1(0), m1(0));

  int K = 3;
  MatrixXd u2 = identity_matrix(K);
  MatrixXd m2(K, K);
  m2 << 1, 0, 0, 0, 1, 0, 0, 0, 1;

  EXPECT_EQ(K, u2.rows());
  EXPECT_EQ(K, u2.cols());
  for (int i = 0; i < m2.size(); i++) {
    EXPECT_EQ(u2(i), m2(i));
  }
}

TEST(MathFunctions, identity_matrix_throw) {
  using stan::math::identity_matrix;
  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(identity_matrix(-1), std::domain_error);
  EXPECT_THROW(identity_matrix(inf), std::domain_error);
  EXPECT_THROW(identity_matrix(nan), std::domain_error);
}
