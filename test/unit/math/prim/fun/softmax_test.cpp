#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathMatrixPrimMat, softmax) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::softmax;

  Matrix<double, Dynamic, 1> x(1);
  x << 0.0;

  Matrix<double, Dynamic, 1> theta = softmax(x);
  EXPECT_EQ(1, theta.size());
  EXPECT_FLOAT_EQ(1.0, theta[0]);

  Matrix<double, Dynamic, 1> x2(2);
  x2 << -1.0, 1.0;
  Matrix<double, Dynamic, 1> theta2 = softmax(x2);
  EXPECT_EQ(2, theta2.size());
  EXPECT_FLOAT_EQ(exp(-1) / (exp(-1) + exp(1)), theta2[0]);
  EXPECT_FLOAT_EQ(exp(1) / (exp(-1) + exp(1)), theta2[1]);

  Matrix<double, Dynamic, 1> x3(3);
  x3 << -1.0, 1.0, 10.0;
  Matrix<double, Dynamic, 1> theta3 = softmax(x3);
  EXPECT_EQ(3, theta3.size());
  EXPECT_FLOAT_EQ(exp(-1) / (exp(-1) + exp(1) + exp(10.0)), theta3[0]);
  EXPECT_FLOAT_EQ(exp(1) / (exp(-1) + exp(1) + exp(10.0)), theta3[1]);
  EXPECT_FLOAT_EQ(exp(10) / (exp(-1) + exp(1) + exp(10.0)), theta3[2]);
}

TEST(MathMatrixPrimMat, softmax_neg_inf) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::softmax;
  constexpr double neg_inf = -std::numeric_limits<double>::infinity();

  // -inf in a vector pins that component to exactly 0; the rest renormalize.
  Matrix<double, Dynamic, 1> v(3);
  v << neg_inf, 1.0, 2.0;
  Matrix<double, Dynamic, 1> theta = softmax(v);
  EXPECT_FLOAT_EQ(0.0, theta[0]);
  EXPECT_FLOAT_EQ(exp(1.0) / (exp(1.0) + exp(2.0)), theta[1]);
  EXPECT_FLOAT_EQ(exp(2.0) / (exp(1.0) + exp(2.0)), theta[2]);
  EXPECT_FLOAT_EQ(1.0, theta.sum());

}

TEST(MathMatrixPrimMat, softmax_row_vector) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::softmax;

  Matrix<double, 1, Dynamic> x(1);
  x << 0.0;
  Matrix<double, 1, Dynamic> theta = softmax(x);
  EXPECT_EQ(1, theta.size());
  EXPECT_FLOAT_EQ(1.0, theta[0]);

  Matrix<double, 1, Dynamic> x2(2);
  x2 << -1.0, 1.0;
  Matrix<double, 1, Dynamic> theta2 = softmax(x2);
  EXPECT_EQ(2, theta2.size());
  EXPECT_FLOAT_EQ(exp(-1) / (exp(-1) + exp(1)), theta2[0]);
  EXPECT_FLOAT_EQ(exp(1) / (exp(-1) + exp(1)), theta2[1]);

  Matrix<double, 1, Dynamic> x3(3);
  x3 << -1.0, 1.0, 10.0;
  Matrix<double, 1, Dynamic> theta3 = softmax(x3);
  EXPECT_EQ(3, theta3.size());
  EXPECT_FLOAT_EQ(exp(-1) / (exp(-1) + exp(1) + exp(10.0)), theta3[0]);
  EXPECT_FLOAT_EQ(exp(1) / (exp(-1) + exp(1) + exp(10.0)), theta3[1]);
  EXPECT_FLOAT_EQ(exp(10) / (exp(-1) + exp(1) + exp(10.0)), theta3[2]);
}
