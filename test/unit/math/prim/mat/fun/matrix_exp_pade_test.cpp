#include <stan/math/prim/mat.hpp>
#include <stan/math/prim/mat/fun/matrix_exp_pade.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>
#include <algorithm>
#include <random>

TEST(MathMatrix, matrix_exp_pade_1x1) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m1(1, 1), m2(1, 1);
  m1 << 0;
  m2 << 1;

  expect_matrix_eq(m2, stan::math::matrix_exp_pade(m1));
}

TEST(MathMatrix, matrix_exp_pade_2x2) {
  // example from Moler & Van Loan, 2003
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m1(2, 2), m2(2, 2);

  m1 << -49, 24, -64, 31;
  m2 << -.735759, .551819, -1.471518, 1.103638;

  expect_matrix_eq(m2, stan::math::matrix_exp_pade(m1));
}

TEST(MathMatrix, matrix_exp_pade_3x3) {
  // example from http:// www.sosmath.com/matrix/expo/expo.html
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m1(3, 3), m2(3, 3);
  m1 << 0, 1, 2, 0, 0, -1, 0, 0, 0;
  m2 << 1, 1, 1.5, 0, 1, -1, 0, 0, 1;

  expect_matrix_eq(m2, stan::math::matrix_exp_pade(m1));
}

TEST(MathMatrix, matrix_exp_pade_3x3_2) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m1(3, 3), m2(3, 3);
  m1 << 89, -66, -18, 20, -14, -4, 360, -270, -73;
  m2 << 245.95891, -182.43047, -49.11821, 93.41549, -67.3433, -18.68310,
      842.54120, -631.90590, -168.14036;

  expect_matrix_eq(m2, stan::math::matrix_exp_pade(m1));
}

TEST(MathMatrix, matrix_exp_100x100) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  int size = 100;
  Matrix<double, Dynamic, Dynamic> S = Eigen::MatrixXd::Identity(size, size),
                                   I = Eigen::MatrixXd::Identity(size, size);
  int col1, col2;
  std::random_device rd;
  std::mt19937 mt(rd());
  for (int i = 0; i < 5 * size; i++) {
    col1 = rd() % size;
    col2 = rd() % size;
    while (col1 == col2)
      col2 = rd() % size;
    S.col(col1) += S.col(col2) * std::pow(-1, rd());
  }
  Matrix<double, Dynamic, Dynamic> S_inv = stan::math::mdivide_right(I, S);
  Matrix<double, 1, Dynamic> diag_elements(size);
  diag_elements.setRandom();
  Matrix<double, 1, Dynamic> exp_diag_elements(size);
  exp_diag_elements = stan::math::exp(diag_elements);

  Matrix<double, Dynamic, Dynamic> A;
  A = S * diag_elements.asDiagonal() * S_inv;
  Matrix<double, Dynamic, Dynamic> exp_A, expm_A;
  exp_A = S * exp_diag_elements.asDiagonal() * S_inv;
  expm_A = stan::math::matrix_exp_pade(A);

  double rel_err
      = 1e-6
        * std::max(exp_A.cwiseAbs().maxCoeff(), expm_A.cwiseAbs().maxCoeff());

  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      EXPECT_NEAR(exp_A(i, j), expm_A(i, j), rel_err);
}
