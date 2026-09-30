#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/math/ad_tape_util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

namespace {
Eigen::MatrixXd scale_multiply_test_matrix(int rows, int cols, double offset) {
  Eigen::MatrixXd x(rows, cols);
  for (int j = 0; j < cols; ++j) {
    for (int i = 0; i < rows; ++i) {
      x(i, j) = 0.3 * (i + 1) - 0.2 * (j + 1) + offset;
    }
  }
  return x;
}
}  // namespace

TEST_F(AgradRev, scale_multiply_values_and_gradients) {
  using stan::math::var;
  const int n = 3, k = 4, m = 2;
  var s = 0.7;
  Eigen::Matrix<var, -1, -1> a = scale_multiply_test_matrix(n, k, 0.1);
  Eigen::Matrix<var, -1, -1> b = scale_multiply_test_matrix(k, m, -0.2);
  const stan::test::tape_delta delta;
  Eigen::Matrix<var, -1, -1> c = stan::math::scale_multiply(s, a, b);
  // one callback plus one vari per result coefficient
  EXPECT_EQ(delta.total(), 1 + n * m);
  EXPECT_MATRIX_NEAR(c.val(), 0.7 * a.val() * b.val(), 1e-12);
  stan::math::sum(c).grad();
  const Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(n, m);
  EXPECT_NEAR(s.adj(), (a.val() * b.val()).sum(), 1e-12);
  EXPECT_MATRIX_NEAR(a.adj(), 0.7 * ones * b.val().transpose(), 1e-12);
  EXPECT_MATRIX_NEAR(b.adj(), 0.7 * a.val().transpose() * ones, 1e-12);
}

TEST_F(AgradRev, scale_multiply_double_scale) {
  using stan::math::var;
  const int n = 2, k = 3, m = 2;
  const Eigen::MatrixXd a = scale_multiply_test_matrix(n, k, 0.1);
  Eigen::Matrix<var, -1, -1> b = scale_multiply_test_matrix(k, m, -0.2);
  Eigen::Matrix<var, -1, -1> c = stan::math::scale_multiply(-1.5, a, b);
  EXPECT_MATRIX_NEAR(c.val(), -1.5 * a * b.val(), 1e-12);
  stan::math::sum(c).grad();
  EXPECT_MATRIX_NEAR(b.adj(),
                     -1.5 * a.transpose() * Eigen::MatrixXd::Ones(n, m), 1e-12);
}
