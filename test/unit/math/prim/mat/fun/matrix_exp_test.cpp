#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <random>

TEST(MathMatrix, matrix_exp_1x1) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m1(1, 1), m2(1, 1);
  m1 << 0;
  m2 << 1;

  expect_matrix_eq(m2, stan::math::matrix_exp(m1));
}

TEST(MathMatrix, matrix_exp_2x2) {
  // example from Moler & Van Loan, 2003
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m1(2, 2), m2(2, 2);
  m1 << -49, 24, -64, 31;
  m2 << -.735759, .551819, -1.471518, 1.103638;

  expect_matrix_eq(m2, stan::math::matrix_exp(m1));
}

TEST(MathMatrix, matrix_exp_2x2_2) {
  // make sure matrix_exp doesn't use matrix_exp_2x2,
  // which would return NaN for this matrix
  // Compare to result from http:// comnuan.com/cmnn01015/
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m(2, 2), exp_m(2, 2);

  m << -0.999984, 0.511211, -0.736924, -0.0826997;

  exp_m << 0.2746483852, 0.2893267425, -0.4170720513, 0.7937977746;

  expect_matrix_eq(exp_m, stan::math::matrix_exp(m));
}

TEST(MathMatrix, matrix_exp_3x3) {
  // example from http:// www.sosmath.com/matrix/expo/expo.html
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m1(3, 3), m2(3, 3);
  m1 << 0, 1, 2, 0, 0, -1, 0, 0, 0;
  m2 << 1, 1, 1.5, 0, 1, -1, 0, 0, 1;

  expect_matrix_eq(m2, stan::math::matrix_exp(m1));
}

TEST(MathMatrix, matrix_exp_3x3_2) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m1(3, 3), m2(3, 3);
  m1 << 89, -66, -18, 20, -14, -4, 360, -270, -73;
  m2 << 245.95891, -182.43047, -49.11821, 93.41549, -67.3433, -18.68310,
      842.54120, -631.90590, -168.14036;

  expect_matrix_eq(m2, stan::math::matrix_exp(m1));
}

TEST(MathMatrix, matrix_exp_100x100) {
  int size = 100;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> S
      = Eigen::MatrixXd::Identity(size, size),
      I = Eigen::MatrixXd::Identity(size, size);
  std::random_device rd;
  std::mt19937 mt(rd());
  int col1, col2;
  for (int i = 0; i < 5 * size; i++) {
    col1 = rd() % size;
    col2 = rd() % size;
    while (col1 == col2)
      col2 = rd() % size;
    S.col(col1) += S.col(col2) * std::pow(-1, rd());
  }
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> S_inv
      = stan::math::mdivide_right(I, S);
  Eigen::Matrix<double, 1, Eigen::Dynamic> diag_elements(size);
  diag_elements.setRandom();
  Eigen::Matrix<double, 1, Eigen::Dynamic> exp_diag_elements
      = stan::math::exp(diag_elements);

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A
      = S * diag_elements.asDiagonal() * S_inv,
      exp_A = S * exp_diag_elements.asDiagonal() * S_inv,
      expm_A = stan::math::matrix_exp(A);

  double rel_err
      = 1e-6
        * std::max(exp_A.cwiseAbs().maxCoeff(), expm_A.cwiseAbs().maxCoeff());

  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      EXPECT_NEAR(exp_A(i, j), expm_A(i, j), rel_err);
}

TEST(MathMatrix, matrix_exp_exceptions) {
  using stan::math::matrix_exp;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m1(0, 0), m2(1, 2);

  EXPECT_THROW(matrix_exp(m1), std::invalid_argument);
  EXPECT_THROW(matrix_exp(m2), std::invalid_argument);
}

TEST(MathMatrix, NOT_A_TEST_matrix_num_err) {
  // Code to showcase how dealing with very small
  // numbers (< 1e-10) can increase the relative
  // error. That is why the conditions for small
  // numbers are laxed (results agree within 1e-10,
  // as oppose to using relative error).

  using stan::math::mdivide_right;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m(2, 2), exp_m(2, 2),
      D(2, 2), A(2, 2), exp_A(2, 2), D_m(2, 2), D_expm(2, 2);
  m << 1e-13, 0, 0, 1e-15;
  D << 1, 2, 3, 4;
  exp_m << exp(1e-13), 0, 0, exp(1e-15);
  D_m = D * m;
  A = mdivide_right(D_m, D);
  D_expm = D * exp_m;

  /*
  std::cout << std::endl;
  std::cout << mdivide_right(D_expm, D) << std::endl << std::endl;
  std::cout << stan::math::matrix_exp(A) << std::endl << std::endl;
  std::cout << stan::math::matrix_exp_pade(A) << std::endl << std::endl;
  */
}
