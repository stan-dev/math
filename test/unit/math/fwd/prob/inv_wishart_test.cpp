#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/special_functions/digamma.hpp>

TEST(ProbDistributionsInvWishart, fvar_double) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::inv_wishart_lpdf;

  Matrix<fvar<double>, Dynamic, Dynamic> Y(3, 3);
  Y << 12.147233, -11.9036079, 1.091046, -11.9036079, 16.7585782, 0.8530256,
      1.091046, 0.8530256, 2.5786609;

  Matrix<fvar<double>, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 7.785215, 3.059788, 1.107166, 3.059788, 10.3515035, -0.1232598,
      1.107166, -0.1232598, 7.7623386;

  double dof = 4.0;
  double log_p = log(2.008407e-08);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      Y(i, j).d_ = 1.0;
      Sigma(i, j).d_ = 1.0;
    }

  EXPECT_NEAR(log_p, stan::math::inv_wishart_lpdf(Y, dof, Sigma).val_, 0.01);
  EXPECT_NEAR(-1.4893348387330674,
              stan::math::inv_wishart_lpdf(Y, dof, Sigma).d_, 0.01);
}

TEST(ProbDistributionsInvWishart, fvar_fvar_double) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::inv_wishart_lpdf;

  Matrix<fvar<fvar<double> >, Dynamic, Dynamic> Y(3, 3);
  Y << 12.147233, -11.9036079, 1.0910458, -11.9036079, 16.7585782, 0.8530256,
      1.0910458, 0.8530256, 2.5786609;

  Matrix<fvar<fvar<double> >, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 7.785215, 3.059788, 1.107166, 3.059788, 10.3515035, -0.1232598,
      1.107166, -0.1232598, 7.7623386;

  double dof = 4.0;
  double log_p = log(2.008407e-08);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      Y(i, j).d_ = 1.0;
      Sigma(i, j).d_ = 1.0;
    }

  EXPECT_NEAR(log_p, stan::math::inv_wishart_lpdf(Y, dof, Sigma).val_.val_,
              0.01);
  EXPECT_NEAR(-1.4893348387330674,
              stan::math::inv_wishart_lpdf(Y, dof, Sigma).d_.val_, 0.01);
}
