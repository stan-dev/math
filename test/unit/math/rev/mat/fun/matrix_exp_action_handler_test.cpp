#include <gtest/gtest.h>
// #include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <test/unit/math/prim/mat/util.hpp>
#include <stan/math/prim/mat/fun/matrix_exp.hpp>
#include <stan/math/prim/mat/fun/matrix_exp_action_handler.hpp>
#include <vector>

TEST(MathMatrix, matrix_exp_action_diag) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  stan::math::matrix_exp_action_handler handler;

  {
    double t = 0.5;
    MatrixXd m1(2, 2);
    VectorXd b(2);
    m1 << 2, 0, 0, 2;
    b << 1, 1;
    auto res = handler.action(m1, b, t);
    EXPECT_NEAR(res(0), M_E, 1.e-8);
    EXPECT_NEAR(res(1), M_E, 1.e-8);
  }

  {
    double t = 1.0;
    MatrixXd m1(2, 2);
    VectorXd b = VectorXd::Random(2);
    m1 << 1, 0, 0, 2;
    auto res = handler.action(m1, b, t);
    EXPECT_NEAR(res(0), b(0) * M_E, 1.e-8);
    EXPECT_NEAR(res(1), b(1) * M_E * M_E, 1.e-8);
  }

  {
    double t = 1.0;
    std::srand(1299);
    MatrixXd m1(2, 2);
    VectorXd b = VectorXd::Random(2);
    m1 << -4.0, 0, 0, -5.0;
    auto res = handler.action(m1, b, t);
    EXPECT_NEAR(res(0), b(0) / (M_E * M_E * M_E * M_E), 1.e-8);
    EXPECT_NEAR(res(1), b(1) / (M_E * M_E * M_E * M_E * M_E), 1.e-8);
  }

  {
    std::srand(999);
    double t = static_cast<double>((std::rand()) / RAND_MAX);
    VectorXd b = VectorXd::Random(5);
    VectorXd d = VectorXd::Random(5);
    MatrixXd m = d.asDiagonal();
    auto res = handler.action(m, b, t);
    for (int i = 0; i < 5; ++i) {
      EXPECT_NEAR(res(i), b(i) * std::exp(t * d(i)), 1.e-8);
    }
  }
}

TEST(MathMatrix, matrix_exp_action_vector) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  stan::math::matrix_exp_action_handler handler;
  std::srand(999);

  for (size_t n = 2; n < 10; ++n) {
    MatrixXd A = MatrixXd::Random(n, n);
    VectorXd b = VectorXd::Random(n);
    VectorXd res = handler.action(A, b);
    MatrixXd expA = stan::math::matrix_exp(A);
    VectorXd expb = expA * b;
    for (size_t i = 0; i < n; ++i) {
      EXPECT_NEAR(res(i), expb(i), 1.e-6);
    }

    int m1, s1, m2, s2;
    const double t1 = 9.9, t2 = 1.0;
    handler.set_approximation_parameter(A, t1, m1, s1);
    A *= t1;
    handler.set_approximation_parameter(A, t2, m2, s2);
    EXPECT_EQ(m1, m2);
    EXPECT_EQ(s1, s2);
  }
}

TEST(MathMatrix, matrix_exp_action_matrix) {
  using Eigen::MatrixXd;
  stan::math::matrix_exp_action_handler handler;
  std::srand(999);

  constexpr int N = 10;
  constexpr int M = 4;
  Eigen::Matrix<double, N, N> A = Eigen::Matrix<double, N, N>::Random();
  Eigen::Matrix<double, N, M> B = Eigen::Matrix<double, N, M>::Random();

  Eigen::Matrix<double, N, M> res = handler.action(A, B);
  MatrixXd Ad(A);
  MatrixXd expa = stan::math::matrix_exp(Ad);
  Eigen::Matrix<double, N, M> expb = expa * B;

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {
      EXPECT_FLOAT_EQ(res(i, j), expb(i, j));
    }
  }
}

TEST(MathMatrix, matrix_exp_action_matrix_transpose) {
  using Eigen::MatrixXd;
  stan::math::matrix_exp_action_handler handler;
  std::srand(1999);

  constexpr int N = 10;
  constexpr int M = 4;
  Eigen::Matrix<double, N, N> A = Eigen::Matrix<double, N, N>::Random();
  Eigen::Matrix<double, N, M> B = Eigen::Matrix<double, N, M>::Random();

  Eigen::Matrix<double, N, M> res = handler.action(A.transpose(), B);
  MatrixXd Ad(A);
  MatrixXd expa = stan::math::matrix_exp(Ad).transpose();
  Eigen::Matrix<double, N, M> expb = expa * B;

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {
      EXPECT_NEAR(res(i, j), expb(i, j), 1.e-6);
    }
  }
}
