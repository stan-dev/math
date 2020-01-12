#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

void test_log_softmax(const Eigen::Matrix<double, Eigen::Dynamic, 1>& theta) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::log_softmax;
  using stan::math::softmax;
  using std::log;

  int theta_size = theta.size();

  Matrix<double, Dynamic, 1> log_softmax_theta = log_softmax(theta);

  Matrix<double, Dynamic, 1> softmax_theta = softmax(theta);

  Matrix<double, Dynamic, 1> log_softmax_theta_expected(theta_size);
  for (int i = 0; i < theta_size; ++i)
    log_softmax_theta_expected(i) = log(softmax_theta(i));

  EXPECT_EQ(log_softmax_theta_expected.size(), log_softmax_theta.size());
  for (int i = 0; i < theta.size(); ++i)
    EXPECT_FLOAT_EQ(log_softmax_theta_expected(i), log_softmax_theta(i));
}

TEST(MathMatrixPrimMat, log_softmax) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::log_softmax;
  using stan::math::softmax;

  // stan::math::vector_d x(1);
  // x << 0.0;
  // test_log_softmax(x);

  stan::math::vector_d x2(2);
  x2 << -1.0, 1.0;
  test_log_softmax(x2);

  // stan::math::vector_d x3(3);
  // x3 << -1.0, 1.0, 10.0;
  // test_log_softmax(x3);
}
TEST(MathMatrixPrimMat, log_softmax_exception) {
  using stan::math::log_softmax;
  stan::math::vector_d v0;  // size == 0

  EXPECT_THROW(log_softmax(v0), std::invalid_argument);
}
