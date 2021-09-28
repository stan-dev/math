#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/distributions.hpp>

TEST(ProbDistributionsWishart, fvar_double) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  Matrix<fvar<double>, Dynamic, Dynamic> Sigma(2, 2);
  Sigma << 1.848220, 1.899623, 1.899623, 12.751941;

  Matrix<fvar<double>, Dynamic, Dynamic> Y(2, 2);
  Y << 2.011108, -11.206611, -11.206611, 112.94139;

  for (int i = 0; i < 4; i++) {
    Sigma(i).d_ = 1.0;
    Y(i).d_ = 1.0;
  }

  unsigned int dof = 3;

  // computed with MCMCpack in R
  double lp = log(8.658e-07);

  EXPECT_NEAR(lp, stan::math::wishart_log(Y, dof, Sigma).val_, 0.01);
  EXPECT_NEAR(-0.76893887, stan::math::wishart_log(Y, dof, Sigma).d_, 0.01);
}

TEST(ProbDistributionsWishart, fvar_fvar_double) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  Matrix<fvar<fvar<double> >, Dynamic, Dynamic> Sigma(2, 2);
  Sigma << 1.848220, 1.899623, 1.899623, 12.751941;

  Matrix<fvar<fvar<double> >, Dynamic, Dynamic> Y(2, 2);
  Y << 2.011108, -11.206611, -11.206611, 112.94139;

  for (int i = 0; i < 4; i++) {
    Sigma(i).d_.val_ = 1.0;
    Y(i).d_.val_ = 1.0;
  }

  unsigned int dof = 3;

  // computed with MCMCpack in R
  double lp = log(8.658e-07);

  EXPECT_NEAR(lp, stan::math::wishart_log(Y, dof, Sigma).val_.val_, 0.01);
  EXPECT_NEAR(-0.76893887, stan::math::wishart_log(Y, dof, Sigma).d_.val_,
              0.01);
}
