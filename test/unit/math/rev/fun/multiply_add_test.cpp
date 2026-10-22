#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

namespace {
Eigen::MatrixXd multiply_add_test_matrix(int rows, int cols, double offset) {
  Eigen::MatrixXd x(rows, cols);
  for (int j = 0; j < cols; ++j) {
    for (int i = 0; i < rows; ++i) {
      x(i, j) = 0.3 * (i + 1) - 0.2 * (j + 1) + offset;
    }
  }
  return x;
}
}  // namespace

TEST_F(AgradRev, multiply_add_values_and_gradients) {
  using stan::math::var;
  const int n = 3, k = 4, m = 2;
  Eigen::Matrix<var, -1, -1> a = multiply_add_test_matrix(n, k, 0.1);
  Eigen::Matrix<var, -1, -1> b = multiply_add_test_matrix(k, m, -0.2);
  Eigen::Matrix<var, -1, -1> c = multiply_add_test_matrix(n, m, 0.4);
  const auto& stack = *stan::math::ChainableStack::instance_;
  const std::size_t start
      = stack.var_stack_.size() + stack.var_nochain_stack_.size();
  Eigen::Matrix<var, -1, -1> d = stan::math::multiply_add(a, b, c);
  // one callback plus one vari per result coefficient
  EXPECT_EQ(stack.var_stack_.size() + stack.var_nochain_stack_.size() - start,
            1 + n * m);
  EXPECT_MATRIX_NEAR(d.val(), a.val() * b.val() + c.val(), 1e-12);
  stan::math::sum(d).grad();
  const Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(n, m);
  EXPECT_MATRIX_NEAR(a.adj(), ones * b.val().transpose(), 1e-12);
  EXPECT_MATRIX_NEAR(b.adj(), a.val().transpose() * ones, 1e-12);
  EXPECT_MATRIX_NEAR(c.adj(), ones, 1e-12);
}

TEST_F(AgradRev, multiply_add_double_addend) {
  using stan::math::var;
  const int n = 2, k = 3, m = 2;
  Eigen::Matrix<var, -1, -1> a = multiply_add_test_matrix(n, k, 0.1);
  const Eigen::MatrixXd b = multiply_add_test_matrix(k, m, -0.2);
  const Eigen::MatrixXd c = multiply_add_test_matrix(n, m, 0.4);
  Eigen::Matrix<var, -1, -1> d = stan::math::multiply_add(a, b, c);
  EXPECT_MATRIX_NEAR(d.val(), a.val() * b + c, 1e-12);
  stan::math::sum(d).grad();
  EXPECT_MATRIX_NEAR(a.adj(), Eigen::MatrixXd::Ones(n, m) * b.transpose(),
                     1e-12);
}
