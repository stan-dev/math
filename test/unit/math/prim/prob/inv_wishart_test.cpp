#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/special_functions/digamma.hpp>

using Eigen::Dynamic;
using Eigen::Matrix;

using stan::math::inv_wishart_log;

TEST(ProbDistributionsInvWishart, LowerTriangular) {
  // Tests if only of the lower triangular portion of
  // outcome and scale matrices are taken
  using Eigen::MatrixXd;
  using stan::math::inv_wishart_log;

  MatrixXd Sigma(4, 4);
  MatrixXd Sigma_sym(4, 4);
  MatrixXd Y(4, 4);
  MatrixXd Y_sym(4, 4);

  Y << 7.988168, -10.955605, -14.47483, 4.395895, -9.555605, 44.750570,
      49.21577, -15.454186, -14.474830, 49.215769, 60.08987, -20.481079,
      4.395895, -18.454186, -21.48108, 7.885833;

  Y_sym << 7.988168, -9.555605, -14.474830, 4.395895, -9.555605, 44.750570,
      49.215769, -18.454186, -14.474830, 49.215769, 60.08987, -21.48108,
      4.395895, -18.454186, -21.48108, 7.885833;

  Sigma << 2.9983662, 0.2898776, -2.650523, 0.1055911, 0.2898776, 11.4803610,
      7.157993, -3.1129955, -2.6505229, 7.1579931, 11.676181, -3.5866852,
      0.1055911, -3.1129955, -3.586685, 1.4482736;

  Sigma_sym << 2.9983662, 0.2898776, -2.6505229, 0.1055911, 0.2898776,
      11.4803610, 7.1579931, -3.1129955, -2.6505229, 7.1579931, 11.676181,
      -3.586685, 0.1055911, -3.1129955, -3.586685, 1.4482736;

  unsigned int dof = 5;

  EXPECT_EQ(inv_wishart_log(Y, dof, Sigma), inv_wishart_log(Y_sym, dof, Sigma));
  EXPECT_EQ(inv_wishart_log(Y, dof, Sigma), inv_wishart_log(Y, dof, Sigma_sym));
  EXPECT_EQ(inv_wishart_log(Y, dof, Sigma),
            inv_wishart_log(Y_sym, dof, Sigma_sym));
}
TEST(ProbDistributionsInvWishart, InvWishart) {
  Matrix<double, Dynamic, Dynamic> Y(3, 3);
  Y << 12.147233, -11.9036079, 1.0910458, -11.903608, 16.7585782, 0.8530256,
      1.091046, 0.8530256, 2.5786609;

  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 7.785215, 3.0597878, 1.1071663, 3.059788, 10.3515035, -0.1232598,
      1.107166, -0.1232598, 7.7623386;

  double dof = 4.0;
  double log_p = log(2.008407e-08);

  EXPECT_NEAR(log_p, stan::math::inv_wishart_log(Y, dof, Sigma), 0.01);
}
TEST(ProbDistributionsInvWishart, Propto) {
  Matrix<double, Dynamic, Dynamic> Y(3, 3);
  Y << 12.147233, -11.9036079, 1.0910458, -11.903608, 16.7585782, 0.8530256,
      1.091046, 0.8530256, 2.5786609;

  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 7.785215, 3.0597878, 1.1071663, 3.059788, 10.3515035, -0.1232598,
      1.107166, -0.1232598, 7.7623386;

  double dof = 4.0;

  EXPECT_FLOAT_EQ(0.0, stan::math::inv_wishart_log<true>(Y, dof, Sigma));
}
TEST(ProbDistributionsInvWishart, Error) {
  Matrix<double, Dynamic, Dynamic> Sigma;
  Matrix<double, Dynamic, Dynamic> Y;
  double nu;

  Sigma.resize(1, 1);
  Y.resize(1, 1);
  Sigma.setIdentity();
  Y.setIdentity();
  nu = 1;
  EXPECT_NO_THROW(inv_wishart_log(Y, nu, Sigma));

  nu = 5;
  Sigma.resize(2, 1);
  EXPECT_THROW(inv_wishart_log(Y, nu, Sigma), std::invalid_argument);

  nu = 5;
  Sigma.resize(2, 2);
  Y.resize(2, 1);
  EXPECT_THROW(inv_wishart_log(Y, nu, Sigma), std::invalid_argument);

  nu = 5;
  Sigma.resize(2, 2);
  Y.resize(3, 3);
  EXPECT_THROW(inv_wishart_log(Y, nu, Sigma), std::invalid_argument);

  Sigma.resize(3, 3);
  Sigma.setIdentity();
  Y.resize(3, 3);
  Y.setIdentity();
  nu = 3;
  EXPECT_NO_THROW(inv_wishart_log(Y, nu, Sigma));
  nu = 2;
  EXPECT_THROW(inv_wishart_log(Y, nu, Sigma), std::domain_error);
}
