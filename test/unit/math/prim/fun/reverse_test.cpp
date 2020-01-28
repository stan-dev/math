#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, reverse_array) {
  using stan::math::reverse;

  std::vector<double> x0;
  EXPECT_EQ(0, reverse(x0).size());

  std::vector<double> x1(5);
  EXPECT_FLOAT_EQ(x1[0], reverse(x1)[0]);

  std::vector<double> x{7, -1.2, 4, -3};
  std::vector<double> rev{-3, 4, -1.2, 7};
  std::vector<double> y = stan::math::reverse(x);
  for (int i = 0; i < x.size(); i++) {
    EXPECT_FLOAT_EQ(y[i], rev[i]);
  }
}

TEST(MathFunctions, reverse_vector) {
  using stan::math::reverse;

  Eigen::VectorXd x0(0);
  EXPECT_EQ(0, reverse(x0).size());

  Eigen::VectorXd x1 = Eigen::VectorXd::Constant(1, 5);
  expect_matrix_eq(x1, reverse(x1));

  Eigen::VectorXd x(4);
  x << 7, -1.2, 4, -3;
  Eigen::VectorXd rev(4);
  rev << -3, 4, -1.2, 7;
  expect_matrix_eq(rev, reverse(x));
}

TEST(MathFunctions, reverse_row_vector) {
  using stan::math::reverse;

  Eigen::RowVectorXd x0(0);
  EXPECT_EQ(0, reverse(x0).size());

  Eigen::RowVectorXd x1 = Eigen::RowVectorXd::Constant(1, 5);
  expect_matrix_eq(x1, reverse(x1));

  Eigen::RowVectorXd x(4);
  x << 7, -1.2, 4, -3;
  Eigen::RowVectorXd rev(4);
  rev << -3, 4, -1.2, 7;
  expect_matrix_eq(rev, reverse(x));
}
