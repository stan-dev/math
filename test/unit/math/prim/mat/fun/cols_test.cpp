#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(primMatFunCols, matrix) {
  using Eigen::MatrixXd;
  using stan::math::cols;
  MatrixXd m(3, 2);
  EXPECT_EQ(2, cols(m));

  MatrixXd m2(3, 1);
  EXPECT_EQ(1, cols(m2));
}
TEST(primMatFunCols, vector) {
  using Eigen::VectorXd;
  using stan::math::cols;
  VectorXd v(2);
  EXPECT_EQ(1, cols(v));

  VectorXd u(1);
  EXPECT_EQ(1, cols(u));
}
TEST(primMatFunCols, rowVector) {
  using Eigen::RowVectorXd;
  using stan::math::cols;
  RowVectorXd v(2);
  EXPECT_EQ(2, cols(v));

  RowVectorXd u(1);
  EXPECT_EQ(1, cols(u));
}
