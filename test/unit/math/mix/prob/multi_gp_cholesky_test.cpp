#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>

TEST_F(AgradRev, ProbDistributionsMultiGPCholesky_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  Matrix<fvar<var>, Dynamic, 1> mu(5, 1);
  mu.setZero();

  Matrix<fvar<var>, Dynamic, Dynamic> y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0, 11.0, 2.0, -5.0, 11.0, 0.0, -2.0, 11.0, 2.0,
      -2.0, -11.0;

  Matrix<fvar<var>, Dynamic, Dynamic> Sigma(5, 5);
  Sigma << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0,
      1.0, 0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;

  Matrix<fvar<var>, Dynamic, 1> w(3, 1);
  w << 1.0, 0.5, 1.5;

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

  Matrix<fvar<var>, Dynamic, Dynamic> L = Sigma.llt().matrixL();

  fvar<var> lp_ref(0);
  for (size_t i = 0; i < 3; i++) {
    Matrix<fvar<var>, Dynamic, 1> cy(y.row(i).transpose());
    Matrix<fvar<var>, Dynamic, Dynamic> cSigma((1.0 / w[i]) * Sigma);
    lp_ref += stan::math::multi_normal_lpdf(cy, mu, cSigma);
  }

  EXPECT_FLOAT_EQ(lp_ref.val_.val(),
                  stan::math::multi_gp_cholesky_lpdf(y, L, w).val_.val());
  EXPECT_FLOAT_EQ(-74.572952,
                  stan::math::multi_gp_cholesky_lpdf(y, L, w).d_.val());

  stan::math::recover_memory();
}

TEST_F(AgradRev, ProbDistributionsMultiGPCholesky_fvar_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  Matrix<fvar<fvar<var> >, Dynamic, 1> mu(5, 1);
  mu.setZero();

  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0, 11.0, 2.0, -5.0, 11.0, 0.0, -2.0, 11.0, 2.0,
      -2.0, -11.0;

  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> Sigma(5, 5);
  Sigma << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0,
      1.0, 0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;

  Matrix<fvar<fvar<var> >, Dynamic, 1> w(3, 1);
  w << 1.0, 0.5, 1.5;

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

  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> L = Sigma.llt().matrixL();

  fvar<fvar<var> > lp_ref(0);
  for (size_t i = 0; i < 3; i++) {
    Matrix<fvar<fvar<var> >, Dynamic, 1> cy(y.row(i).transpose());
    Matrix<fvar<fvar<var> >, Dynamic, Dynamic> cSigma((1.0 / w[i]) * Sigma);
    lp_ref += stan::math::multi_normal_lpdf(cy, mu, cSigma);
  }

  EXPECT_FLOAT_EQ(lp_ref.val_.val_.val(),
                  stan::math::multi_gp_cholesky_lpdf(y, L, w).val_.val_.val());
  EXPECT_FLOAT_EQ(-74.572952,
                  stan::math::multi_gp_cholesky_lpdf(y, L, w).d_.val_.val());

  stan::math::recover_memory();
}
