#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>

TEST(ProbDistributionsMultiNormal, fvar_double) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using std::vector;

  Matrix<fvar<double>, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<fvar<double>, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<fvar<double>, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  for (int i = 0; i < 3; i++) {
    y(i).d_ = 1.0;
    mu(i).d_ = 1.0;
    for (int j = 0; j < 3; j++)
      Sigma(i, j).d_ = 1.0;
  }
  EXPECT_FLOAT_EQ(-11.73908, stan::math::multi_normal_lpdf(y, mu, Sigma).val_);
  EXPECT_FLOAT_EQ(0.54899865, stan::math::multi_normal_lpdf(y, mu, Sigma).d_);
}

TEST(ProbDistributionsMultiNormal, fvar_fvar_double) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using std::vector;

  Matrix<fvar<fvar<double> >, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<fvar<fvar<double> >, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<fvar<fvar<double> >, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  for (int i = 0; i < 3; i++) {
    y(i).d_ = 1.0;
    mu(i).d_ = 1.0;
    for (int j = 0; j < 3; j++)
      Sigma(i, j).d_ = 1.0;
  }
  EXPECT_FLOAT_EQ(-11.73908,
                  stan::math::multi_normal_lpdf(y, mu, Sigma).val_.val_);
  EXPECT_FLOAT_EQ(0.54899865,
                  stan::math::multi_normal_lpdf(y, mu, Sigma).d_.val_);
}
