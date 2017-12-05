#include <stan/math/fwd/mat.hpp>
#include <gtest/gtest.h>

using Eigen::Dynamic;
using Eigen::Matrix;

TEST(ProbDistributionsMultiGP, fvar_double) {
  using stan::math::fvar;
  Matrix<fvar<double>, Dynamic, 1> mu(5, 1);
  mu.setZero();

  Matrix<fvar<double>, Dynamic, Dynamic> y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0,
       11.0, 2.0, -5.0, 11.0, 0.0,
      -2.0, 11.0, 2.0, -2.0, -11.0;

  Matrix<fvar<double>, Dynamic, Dynamic> Sigma(5, 5);
  Sigma << 9.0, -3.0, 0.0,  0.0, 0.0,
          -3.0,  4.0, 0.0,  0.0, 0.0,
           0.0,  0.0, 5.0,  1.0, 0.0,
           0.0,  0.0, 1.0, 10.0, 0.0,
           0.0,  0.0, 0.0,  0.0, 2.0;

  Matrix<fvar<double>, Dynamic, 1> w(3, 1);
  w << 1.0,
       0.5,
       1.5;

  for (int i = 0; i < 5; i++) {
    mu(i).d_ = 1.0;
    if (i < 3)
      w(i).d_ = 1.0;
    for (int j = 0; j < 5; j++) {
      Sigma(i, j).d_ = 1.0;
      if (i < 3)
        y(i, j).d_ = 1.0;
    }
  }

  fvar<double> lp_ref(0);
  for (size_t i = 0; i < 3; i++) {
    Matrix<fvar<double>, Dynamic, 1> cy(y.row(i).transpose());
    Matrix<fvar<double>, Dynamic, Dynamic> cSigma((1.0/w[i])*Sigma);
    lp_ref += stan::math::multi_normal_log(cy, mu, cSigma);
  }

  EXPECT_FLOAT_EQ(lp_ref.val_, stan::math::multi_gp_log(y, Sigma, w).val_);
  EXPECT_FLOAT_EQ(-74.572952, stan::math::multi_gp_log(y, Sigma, w).d_);
}

TEST(ProbDistributionsMultiGP, fvar_fvar_double) {
  using stan::math::fvar;
  Matrix<fvar<fvar<double> >, Dynamic, 1> mu(5, 1);
  mu.setZero();

  Matrix<fvar<fvar<double> >, Dynamic, Dynamic> y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0,
       11.0, 2.0, -5.0, 11.0, 0.0,
      -2.0, 11.0, 2.0, -2.0, -11.0;

  Matrix<fvar<fvar<double> >, Dynamic, Dynamic> Sigma(5, 5);
  Sigma << 9.0, -3.0, 0.0,  0.0, 0.0,
          -3.0,  4.0, 0.0,  0.0, 0.0,
           0.0,  0.0, 5.0,  1.0, 0.0,
           0.0,  0.0, 1.0, 10.0, 0.0,
           0.0,  0.0, 0.0,  0.0, 2.0;

  Matrix<fvar<fvar<double> >, Dynamic, 1> w(3, 1);
  w << 1.0,
       0.5,
       1.5;

  for (int i = 0; i < 5; i++) {
    mu(i).d_.val_ = 1.0;
    if (i < 3)
      w(i).d_.val_ = 1.0;
    for (int j = 0; j < 5; j++) {
      Sigma(i, j).d_.val_ = 1.0;
      if (i < 3)
        y(i, j).d_.val_ = 1.0;
    }
  }

  fvar<fvar<double> > lp_ref(0);
  for (size_t i = 0; i < 3; i++) {
    Matrix<fvar<fvar<double> >, Dynamic, 1> cy(y.row(i).transpose());
    Matrix<fvar<fvar<double> >, Dynamic, Dynamic> cSigma((1.0/w[i])*Sigma);
    lp_ref += stan::math::multi_normal_log(cy, mu, cSigma);
  }

  EXPECT_FLOAT_EQ(lp_ref.val_.val_,
                  stan::math::multi_gp_log(y, Sigma, w).val_.val_);
  EXPECT_FLOAT_EQ(-74.572952, stan::math::multi_gp_log(y, Sigma, w).d_.val_);
}
