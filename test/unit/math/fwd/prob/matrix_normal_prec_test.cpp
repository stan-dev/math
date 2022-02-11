#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>

TEST(ProbDistributionsMatrixNormal, fvar_double) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;

  Matrix<fvar<double>, Dynamic, Dynamic> y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0, 11.0, 2.0, -5.0, 11.0, 0.0, -2.0, 11.0, 2.0,
      -2.0, -11.0;

  Matrix<fvar<double>, Dynamic, Dynamic> mu(3, 5);
  mu.setZero();

  Matrix<fvar<double>, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 1.0, 0.5, 0.1, 0.5, 1.0, 0.2, 0.1, 0.2, 1.0;

  Matrix<fvar<double>, Dynamic, Dynamic> D(5, 5);
  D << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 1.0,
      0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;

  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 5; j++) {
      D(i, j).d_ = 1.0;
      if (i < 3) {
        mu(i, j).d_ = 1.0;
        y(i, j).d_ = 1.0;
        if (j < 3)
          Sigma(i, j).d_ = 1.0;
      }
    }

  fvar<double> lp_ref = stan::math::matrix_normal_prec_lpdf(y, mu, Sigma, D);
  EXPECT_FLOAT_EQ(-2132.07482, lp_ref.val_);
  EXPECT_FLOAT_EQ(-2075.1274, lp_ref.d_);
}

TEST(ProbDistributionsMatrixNormal, fvar_fvar_double) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;

  Matrix<fvar<fvar<double>>, Dynamic, Dynamic> y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0, 11.0, 2.0, -5.0, 11.0, 0.0, -2.0, 11.0, 2.0,
      -2.0, -11.0;

  Matrix<fvar<fvar<double>>, Dynamic, Dynamic> mu(3, 5);
  mu.setZero();

  Matrix<fvar<fvar<double>>, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 1.0, 0.5, 0.1, 0.5, 1.0, 0.2, 0.1, 0.2, 1.0;

  Matrix<fvar<fvar<double>>, Dynamic, Dynamic> D(5, 5);
  D << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 1.0,
      0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;

  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 5; j++) {
      D(i, j).d_.val_ = 1.0;
      if (i < 3) {
        mu(i, j).d_.val_ = 1.0;
        y(i, j).d_.val_ = 1.0;
        if (j < 3)
          Sigma(i, j).d_.val_ = 1.0;
      }
    }

  fvar<fvar<double>> lp_ref
      = stan::math::matrix_normal_prec_lpdf(y, mu, Sigma, D);
  EXPECT_FLOAT_EQ(-2132.07482, lp_ref.val_.val_);
  EXPECT_FLOAT_EQ(-2075.1274, lp_ref.d_.val_);
}
