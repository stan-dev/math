#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

void test_log_softmax(const Eigen::VectorXd& x) {
  using Eigen::VectorXd;
  using stan::math::log_softmax;
  using stan::math::softmax;
  using stan::math::log;
  VectorXd y = log_softmax(x);
  VectorXd y_expected = log(softmax(x));
  EXPECT_EQ(y_expected.size(), y.size());
  for (int i = 0; i < y_expected.size(); ++i)
    EXPECT_FLOAT_EQ(y_expected(i), y(i));
}

TEST(MathMatrixPrimMat, log_softmax) {
  using Eigen::VectorXd;
  using stan::math::log_softmax;
  using stan::math::softmax;
  using stan::math::log;

  VectorXd x1(1);
  x1 << 0.0;
  test_log_softmax(x1);

  VectorXd x2(2);
  x2 << -1, 1;
  test_log_softmax(x2);

  VectorXd x4(4);
  x4 << -1.0, 1.0, -1.0, 1.0;
  test_log_softmax(x4);
}

TEST(MathMatrixPrimMat, log_softmax_exception) {
  using stan::math::log_softmax;
  Eigen::VectorXd v0(0);
  EXPECT_THROW(log_softmax(v0), std::invalid_argument);
}
