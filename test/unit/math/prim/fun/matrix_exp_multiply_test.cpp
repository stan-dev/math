#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/util.hpp>
#include <gtest/gtest.h>
#include <vector>

template <int N, int M>
inline void test_matrix_exp_multiply() {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::matrix_exp_multiply;

  std::srand(1999);

  Eigen::MatrixXd A = Eigen::MatrixXd::Random(N, N);
  Eigen::Matrix<double, -1, M> B = Eigen::MatrixXd::Random(N, M);
  Eigen::Matrix<double, Dynamic, Dynamic> A0 = A;

  // brute force
  Eigen::Matrix<double, N, M> expAB
      = stan::math::multiply(stan::math::matrix_exp(A0), B);

  // matrix_exp_multiply
  Eigen::Matrix<double, N, M> res = matrix_exp_multiply(A, B);

  EXPECT_MATRIX_FLOAT_EQ(expAB, res);
}

TEST(MathMatrixPrimMat, matrix_exp_multiply_power_1_norm_fun) {
  stan::math::matrix_exp_action_handler maexp;
  Eigen::MatrixXd x(2, 2);
  x << 1, -2, -3, 4;
  EXPECT_FLOAT_EQ(maexp.mat_power_1_norm(x, 2), 32.0);
  EXPECT_FLOAT_EQ(maexp.mat_power_1_norm(x, 3), 172.0);

  Eigen::MatrixXd y(3, 3);
  y << 1, 2, 3, 4, 5, 6, 7, 1, 2;
  EXPECT_FLOAT_EQ(maexp.mat_power_1_norm(y, 3), 1163.0);
  EXPECT_FLOAT_EQ(maexp.mat_power_1_norm(y, 12), 8.3805595e+11);
}

TEST(MathMatrixPrimMat, matrix_exp_multiply_approx_order) {
  stan::math::matrix_exp_action_handler maexp;
  Eigen::MatrixXd x(2, 2);
  x << 1, -2, -3, 4;
  Eigen::VectorXd b(2);
  b << 1.5, 2.5;
  int m, s;
  maexp.set_approx_order(x, b, 1.0, m, s);
  EXPECT_EQ(m, 40);
  EXPECT_EQ(s, 1);
}

TEST(MathMatrixPrimMat, matrix_exp_multiply) {
  test_matrix_exp_multiply<1, 1>();
  test_matrix_exp_multiply<1, 5>();
  test_matrix_exp_multiply<5, 1>();
  test_matrix_exp_multiply<5, 5>();
  test_matrix_exp_multiply<20, 2>();
}

TEST(MathMatrixPrimMat, matrix_exp_multiply_issue_2529) {
  // issue #2529 https://github.com/stan-dev/math/issues/2529
  Eigen::MatrixXd a(3, 3);
  a << -1.2, 1.2, 0, 0, -0.5, 0.5, 0, 0, -0.3;
  Eigen::VectorXd b(3);
  b << 1.0, 1.0, 1.0;
  for (auto i = 15; i < 40; ++i) {
    double t = i;
    Eigen::MatrixXd m1 = stan::math::matrix_exp_multiply(a * t, b);
    Eigen::MatrixXd m2 = stan::math::matrix_exp(a * t) * b;
    EXPECT_MATRIX_FLOAT_EQ(m1, m2);
  }
}

TEST(MathMatrixPrimMat, matrix_exp_multiply_poisson_5) {
  // Block tridiag mat by discretizing Poisson's Eq. on a 5x5 grid
  int n = 5;
  int nd = n * n;
  Eigen::MatrixXd m = Eigen::MatrixXd::Zero(nd, nd);
  Eigen::VectorXd b(nd);
  for (auto i = 0; i < nd; ++i) {
    b(i) = -1.0 + i * 2.0 / (nd - 1);
    m(i, i) = 4.0;
    if (i + 1 < nd) {
      m(i, i + 1) = -1.0;
      m(i + 1, i) = -1.0;
    }
    if (i + 5 < nd) {
      m(i, std::min(i + 5, nd - 1)) = -1.0;
      m(std::min(i + 5, nd - 1), i) = -1.0;
    }
  }

  double t = 1.0;
  Eigen::MatrixXd p1 = stan::math::matrix_exp_multiply(m * t, b);
  Eigen::MatrixXd p2 = stan::math::matrix_exp(m * t) * b;
  EXPECT_MATRIX_NEAR(p1, p2, 1.e-12);
}

TEST(MathMatrixPrimMat, matrix_exp_multiply_exception) {
  using stan::math::matrix_exp_multiply;
  {  // multiplicable
    Eigen::MatrixXd A(0, 0);
    Eigen::MatrixXd B = Eigen::MatrixXd::Random(1, 2);
    EXPECT_THROW(matrix_exp_multiply(A, B), std::invalid_argument);
    EXPECT_THROW(matrix_exp_multiply(B, A), std::invalid_argument);
  }

  {  // multiplicable
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(2, 2);
    Eigen::MatrixXd B = Eigen::MatrixXd::Random(3, 2);
    EXPECT_THROW(matrix_exp_multiply(A, B), std::invalid_argument);
  }

  {  // square
    Eigen::MatrixXd A(0, 1);
    Eigen::MatrixXd B(1, 2);
    EXPECT_THROW(matrix_exp_multiply(A, B), std::invalid_argument);
  }

  {  // square
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(2, 3);
    Eigen::MatrixXd B = Eigen::MatrixXd::Random(3, 2);
    EXPECT_THROW(matrix_exp_multiply(A, B), std::invalid_argument);
  }
}
