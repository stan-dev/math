#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/special_functions/digamma.hpp>

TEST(ProbDistributionsInvWishart, inv_wishart_symmetry) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  using Eigen::MatrixXd;
  using stan::math::inv_wishart_lpdf;

  MatrixXd Sigma(4, 4);
  MatrixXd Sigma_non_sym(4, 4);
  MatrixXd Y(4, 4);
  MatrixXd Y_non_sym(4, 4);

  Y << 7.988168, -9.555605, -14.47483, 4.395895, -9.555605, 44.750570,
      49.215769, -15.454186, -14.474830, 49.215769, 60.08987, -20.48108,
      4.395895, -15.454186, -20.48108, 7.885833;

  Sigma << 2.9983662, 0.2898776, -2.650523, 0.1055911, 0.2898776, 11.4803610,
      7.1579931, -3.1129955, -2.650523, 7.1579931, 11.676181, -3.586685,
      0.1055911, -3.1129955, -3.586685, 1.4482736;

  Y_non_sym << 7.988168, 100.0, -14.47483, 4.395895, -9.555605, 44.750570,
      49.21577, -15.454186, -14.474830, 49.215769, 60.08987, -20.481079,
      4.395895, -18.454186, -21.48108, 7.885833;

  Sigma_non_sym << 2.9983662, 100.0, -2.650523, 0.1055911, 0.2898776,
      11.4803610, 7.1579931, -3.1129955, -2.650523, 7.1579931, 11.676181,
      -3.586685, 0.1055911, -3.1129955, -3.586685, 1.4482736;

  unsigned int dof = 5;

  EXPECT_NO_THROW(inv_wishart_lpdf(Y, dof, Sigma));
  EXPECT_THROW(inv_wishart_lpdf(Y_non_sym, dof, Sigma), std::domain_error);
  EXPECT_THROW(inv_wishart_lpdf(Y, dof, Sigma_non_sym), std::domain_error);
}

TEST(ProbDistributionsInvWishart, inv_wishart_pos_def) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  using Eigen::MatrixXd;
  using stan::math::inv_wishart_lpdf;

  MatrixXd Sigma(2, 2);
  MatrixXd Sigma_non_pos_def(2, 2);
  MatrixXd Y(2, 2);
  MatrixXd Y_non_pos_def(2, 2);

  Sigma << 1, 0, 0, 1;
  Sigma_non_pos_def << -1, 0, 0, 1;
  Y << 1, 0, 0, 1;
  Y_non_pos_def << -1, 0, 0, 1;

  unsigned int dof = 5;

  EXPECT_NO_THROW(inv_wishart_lpdf(Y, dof, Sigma));
  EXPECT_THROW(inv_wishart_lpdf(Y_non_pos_def, dof, Sigma), std::domain_error);
  EXPECT_THROW(inv_wishart_lpdf(Y, dof, Sigma_non_pos_def), std::domain_error);
}

TEST(ProbDistributionsInvWishart, InvWishart) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  using stan::math::inv_wishart_lpdf;
  Matrix<double, Dynamic, Dynamic> Y(3, 3);
  Y << 12.147233, -11.903608, 1.0910458, -11.903608, 16.7585782, 0.8530256,
      1.0910458, 0.8530256, 2.5786609;

  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 7.785215, 3.059788, 1.107166, 3.059788, 10.3515035, -0.1232598,
      1.107166, -0.1232598, 7.7623386;

  double dof = 4.0;
  double log_p = log(2.008407e-08);

  EXPECT_NEAR(log_p, stan::math::inv_wishart_lpdf(Y, dof, Sigma), 0.01);
}
TEST(ProbDistributionsInvWishart, Propto) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  using stan::math::inv_wishart_lpdf;
  Matrix<double, Dynamic, Dynamic> Y(3, 3);
  Y << 12.147233, -11.903608, 1.091046, -11.903608, 16.7585782, 0.8530256,
      1.091046, 0.8530256, 2.5786609;

  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 7.785215, 3.059788, 1.107166, 3.059788, 10.3515035, -0.1232598,
      1.107166, -0.1232598, 7.7623386;

  double dof = 4.0;

  EXPECT_FLOAT_EQ(0.0, stan::math::inv_wishart_lpdf<true>(Y, dof, Sigma));
}
TEST(ProbDistributionsInvWishart, Error) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  using stan::math::inv_wishart_lpdf;
  Matrix<double, Dynamic, Dynamic> Sigma;
  Matrix<double, Dynamic, Dynamic> Y;
  double nu;

  Sigma.resize(1, 1);
  Y.resize(1, 1);
  Sigma.setIdentity();
  Y.setIdentity();
  nu = 1;
  EXPECT_NO_THROW(inv_wishart_lpdf(Y, nu, Sigma));

  nu = 5;
  Sigma.resize(2, 1);
  EXPECT_THROW(inv_wishart_lpdf(Y, nu, Sigma), std::invalid_argument);

  nu = 5;
  Sigma.resize(2, 2);
  Y.resize(2, 1);
  EXPECT_THROW(inv_wishart_lpdf(Y, nu, Sigma), std::invalid_argument);

  nu = 5;
  Sigma.resize(2, 2);
  Y.resize(3, 3);
  EXPECT_THROW(inv_wishart_lpdf(Y, nu, Sigma), std::invalid_argument);

  Sigma.resize(3, 3);
  Sigma.setIdentity();
  Y.resize(3, 3);
  Y.setIdentity();
  nu = 3;
  EXPECT_NO_THROW(inv_wishart_lpdf(Y, nu, Sigma));
  nu = 2;
  EXPECT_THROW(inv_wishart_lpdf(Y, nu, Sigma), std::domain_error);
}
