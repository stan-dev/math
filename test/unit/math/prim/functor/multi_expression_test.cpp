#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <tuple>

TEST(MathFunctions, multi_expression) {
  Eigen::MatrixXd a(2, 2);
  a << 1, 2, 3, 4;
  Eigen::MatrixXd b(2, 2);
  b << 1, 2, 3, 5;
  Eigen::MatrixXd c, d;
  stan::math::eigen_results(c, d) = stan::math::eigen_expressions(a * 2, b + a);

  Eigen::MatrixXd c_res = a * 2;
  Eigen::MatrixXd d_res = b + a;
  EXPECT_MATRIX_EQ(c, c_res);
  EXPECT_MATRIX_EQ(d, d_res);
}
