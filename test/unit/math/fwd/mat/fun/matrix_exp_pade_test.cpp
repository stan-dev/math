#include <stan/math/fwd/mat.hpp>
#include <stan/math/prim/mat/fun/matrix_exp_pade.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <random>

using stan::math::fvar;
using stan::math::matrix_exp_pade;
using stan::math::matrix_fd;

TEST(MathMatrix, matrix_exp_pade_1x1) {
  using stan::math::fvar;

  fvar<double> a;
  a.val_ = 0.0;
  a.d_ = 1.0;

  stan::math::matrix_fd input(1, 1), output;
  input << a;
  output = matrix_exp_pade(input);

  EXPECT_EQ(1.0, output(0, 0).val_);
  EXPECT_EQ(1.0, output(0, 0).d_);

  fvar<double> b;
  b.val_ = 1.0;
  b.d_ = 2.0;

  input << b;
  output = stan::math::matrix_exp_pade(input);

  EXPECT_EQ(exp(1.0), output(0, 0).val_);
  EXPECT_EQ(b.d_ * output(0, 0).val_, output(0, 0).d_);
}

TEST(MathMatrix, matrix_exp_pade_2x2) {
  // example from Moler & Van Loan, 2003, section 3
  stan::math::fvar<double> a, b;
  a.val_ = -1.0;
  b.val_ = -17.0;
  a.d_ = 1.0;
  b.d_ = 1.0;

  stan::math::matrix_fd input(2, 2);
  input << -2 * a + 3 * b, 1.5 * a - 1.5 * b, -4 * a + 4 * b, 3 * a - 2 * b;

  stan::math::matrix_fd output;
  output = stan::math::matrix_exp_pade(input);

  EXPECT_FLOAT_EQ(-0.735759, output(0, 0).val_);
  EXPECT_FLOAT_EQ(0.551819, output(0, 1).val_);
  EXPECT_FLOAT_EQ(-1.471518, output(1, 0).val_);
  EXPECT_FLOAT_EQ(1.103638, output(1, 1).val_);

  // note: in this particular example, derivatives
  // is the same as the value, due to the way the
  // input matrix was constructed.
  EXPECT_FLOAT_EQ(-0.735759, output(0, 0).d_);
  EXPECT_FLOAT_EQ(0.551819, output(0, 1).d_);
  EXPECT_FLOAT_EQ(-1.471518, output(1, 0).d_);
  EXPECT_FLOAT_EQ(1.103638, output(1, 1).d_);
}

TEST(MathMatrix, matrix_exp_pade_3x3) {
  using stan::math::fvar;
  using stan::math::matrix_fd;

  stan::math::fvar<double> a, b, c;
  a.val_ = -1.0;
  b.val_ = 2.0;
  c.val_ = 1.0;
  a.d_ = 1.0;
  b.d_ = 1.0;
  c.d_ = 1.0;

  stan::math::matrix_fd input(3, 3);
  input << -24 * a + 40 * b - 15 * c, 18 * a - 30 * b + 12 * c,
      5 * a - 8 * b + 3 * c, 20 * b - 20 * c, -15 * b + 16 * c, -4 * b + 4 * c,
      -120 * a + 120 * b, 90 * a - 90 * b, 25 * a - 24 * b;

  stan::math::matrix_fd output;
  output = matrix_exp_pade(input);

  EXPECT_FLOAT_EQ(245.95891, output(0, 0).val_);
  EXPECT_FLOAT_EQ(-182.43047, output(0, 1).val_);
  EXPECT_FLOAT_EQ(-49.11821, output(0, 2).val_);
  EXPECT_FLOAT_EQ(93.41549, output(1, 0).val_);
  EXPECT_FLOAT_EQ(-67.3433, output(1, 1).val_);
  EXPECT_FLOAT_EQ(-18.68310, output(1, 2).val_);
  EXPECT_FLOAT_EQ(842.54120, output(2, 0).val_);
  EXPECT_FLOAT_EQ(-631.90590, output(2, 1).val_);
  EXPECT_FLOAT_EQ(-168.14036, output(2, 2).val_);

  // Note: because of the way the matrix was constructed,
  // derivative should be the same as value (same case
  // as in the previous test).

  EXPECT_FLOAT_EQ(245.95891, output(0, 0).d_);
  EXPECT_FLOAT_EQ(-182.43047, output(0, 1).d_);
  EXPECT_FLOAT_EQ(-49.11821, output(0, 2).d_);
  EXPECT_FLOAT_EQ(93.41549, output(1, 0).d_);
  EXPECT_FLOAT_EQ(-67.3433, output(1, 1).d_);
  EXPECT_FLOAT_EQ(-18.68310, output(1, 2).d_);
  EXPECT_FLOAT_EQ(842.54120, output(2, 0).d_);
  EXPECT_FLOAT_EQ(-631.90590, output(2, 1).d_);
  EXPECT_FLOAT_EQ(-168.14036, output(2, 2).d_);
}

TEST(MathMatrix, matrix_exp_100x100) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::matrix_fd;

  int size = 100;
  std::random_device rd;
  std::mt19937 mt(rd());
  Matrix<double, Dynamic, Dynamic> S = Eigen::MatrixXd::Identity(size, size),
                                   I = Eigen::MatrixXd::Identity(size, size);
  int col1, col2;
  for (int i = 0; i < 5 * size; i++) {
    col1 = rd() % size;
    col2 = rd() % size;
    while (col1 == col2)
      col2 = rd() % size;
    S.col(col1) += S.col(col2) * std::pow(-1, rd());
  }
  Matrix<double, Dynamic, Dynamic> S_inv = stan::math::mdivide_right(I, S);
  Matrix<double, 1, Dynamic> diag_elements_d(size);
  diag_elements_d.setRandom();
  matrix_fd diag_elements = diag_elements_d.cast<fvar<double> >();
  for (int i = 0; i < size; i++)
    diag_elements(i).d_ = 1.0;

  Matrix<double, Dynamic, Dynamic> exp_diag_elements_d
      = stan::math::exp(diag_elements_d),
      exp_A = S * exp_diag_elements_d.asDiagonal() * S_inv;

  matrix_fd A = S.cast<fvar<double> >() * diag_elements.asDiagonal()
                * S_inv.cast<fvar<double> >(),
            expm_A = stan::math::matrix_exp_pade(A);

  // Note: because of the way of matrix was constructed,
  // derivative should be the same as value (same case
  // as in the previous test).
  double rel_err = 1e-6
                   * std::max(exp_A.cwiseAbs().maxCoeff(),
                              expm_A.cwiseAbs().maxCoeff().val_),
         rel_err_d = 1e-6
                     * std::max(exp_A.cwiseAbs().maxCoeff(),
                                expm_A.cwiseAbs().maxCoeff().d_);

  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      EXPECT_NEAR(exp_A(i, j), expm_A(i, j).val_, rel_err);
      EXPECT_NEAR(exp_A(i, j), expm_A(i, j).d_, rel_err_d);
    }
  }
}
