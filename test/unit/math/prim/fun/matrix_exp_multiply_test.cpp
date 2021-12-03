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

TEST(MathMatrixPrimMat, matrix_exp_multiply_matlab) {
  // https://blogs.mathworks.com/cleve/2012/07/23/a-balancing-act-for-the-matrix-exponential/#7401fe8a-5a7d-40df-92d7-1ae34f45adf2
  // // NOLINT
  {
    Eigen::MatrixXd m(3, 3);
    m << 0, 1e-08, 0, -20066666666.6667, -3, 20000000000, 66.6666666666667, 0,
        -66.6666666666667;
    Eigen::MatrixXd expm(3, 3);
    expm << 0.446849468283175, 1.54044157383952e-09, 0.462811453558774,
        -5743067.77947947, -0.0152830038686819, -4526542.71278401,
        0.447722977849494, 1.54270484519591e-09, 0.463480648837651;
    std::vector<double> r{-1, 0, 1};
    Eigen::VectorXd p1(3), p2(3);
    for (auto i : r) {
      for (auto j : r) {
        for (auto k : r) {
          Eigen::VectorXd b(3);
          b << i, j, k;
          p1 = stan::math::matrix_exp_multiply(m, b);
          p2 = expm * b;
          EXPECT_MATRIX_NEAR(p1, p2, 5.e-6);
        }
      }
    }
  }

  // https://discourse.mc-stan.org/t/documentation-for-the-action-of-the-matrix-exponential/25498/8?u=yizhang
  // // NOLINT
  {
    Eigen::MatrixXd m(2, 2);
    m << -49, 24, -64, 31;
    Eigen::MatrixXd expm(2, 2);
    expm << -0.735758758144755, 0.551819099658099, -1.47151759908826,
        1.10363824071558;
    std::vector<double> r{0, 1, 10, 100, 1000, 10000};
    Eigen::VectorXd p1(2), p2(2);
    for (auto i : r) {
      for (auto j : r) {
        Eigen::VectorXd b(2);
        b << i, j;
        p1 = stan::math::matrix_exp_multiply(m, b);
        p2 = expm * b;
        EXPECT_MATRIX_NEAR(p1, p2, 1.e-8);
      }
    }
  }
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
