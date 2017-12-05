#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>

using Eigen::Dynamic;
using Eigen::Matrix;
using std::vector;

TEST(ProbDistributionsMultiNormalPrec, fvar_var) {
  using stan::math::fvar;
  using stan::math::var;
  Matrix<fvar<var>, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<fvar<var>, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<fvar<var>, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0,
    -3.0,  4.0, 0.0,
    0.0, 0.0, 5.0;

  Sigma << 9.0, -3.0, 0.0,
    -3.0,  4.0, 0.0,
    0.0, 0.0, 5.0;
  for (int i = 0; i < 3; i++) {
    y(i).d_ = 1.0;
    mu(i).d_ = 1.0;
    for (int j = 0; j < 3; j++)
      Sigma(i, j).d_ = 1.0;
  }

  Matrix<fvar<var>, Dynamic, Dynamic> L = Sigma.inverse();
  EXPECT_FLOAT_EQ(-11.73908,
                  stan::math::multi_normal_prec_log(y, mu, L).val_.val());
  EXPECT_FLOAT_EQ(0.54899865,
                  stan::math::multi_normal_prec_log(y, mu, L).d_.val());
}

TEST(ProbDistributionsMultiNormalPrec, fvar_fvar_var) {
  using stan::math::fvar;
  using stan::math::var;
  Matrix<fvar<fvar<var> >, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<fvar<fvar<var> >, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0,
    -3.0,  4.0, 0.0,
    0.0, 0.0, 5.0;

  Sigma << 9.0, -3.0, 0.0,
    -3.0,  4.0, 0.0,
    0.0, 0.0, 5.0;
  for (int i = 0; i < 3; i++) {
    y(i).d_.val_ = 1.0;
    mu(i).d_.val_ = 1.0;
    for (int j = 0; j < 3; j++)
      Sigma(i, j).d_.val_ = 1.0;
  }

  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> L = Sigma.inverse();
  EXPECT_FLOAT_EQ(-11.73908,
                  stan::math::multi_normal_prec_log(y, mu, L).val_.val_.val());
  EXPECT_FLOAT_EQ(0.54899865,
                  stan::math::multi_normal_prec_log(y, mu, L).d_.val_.val());
}
