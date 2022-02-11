#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/distributions.hpp>

TEST(ProbDistributionsWishart, wishart_symmetry) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  using Eigen::MatrixXd;
  using stan::math::wishart_log;

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

  EXPECT_NO_THROW(wishart_log(Y, dof, Sigma));
  EXPECT_THROW(wishart_log(Y_non_sym, dof, Sigma), std::domain_error);
  EXPECT_THROW(wishart_log(Y, dof, Sigma_non_sym), std::domain_error);
}

TEST(ProbDistributionsWishart, wishart_pos_def) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  using Eigen::MatrixXd;
  using stan::math::wishart_log;

  MatrixXd Sigma(2, 2);
  MatrixXd Sigma_non_pos_def(2, 2);
  MatrixXd Y(2, 2);
  MatrixXd Y_non_pos_def(2, 2);

  Sigma << 1, 0, 0, 1;
  Sigma_non_pos_def << -1, 0, 0, 1;
  Y << 1, 0, 0, 1;
  Y_non_pos_def << -1, 0, 0, 1;

  unsigned int dof = 5;

  EXPECT_NO_THROW(wishart_log(Y, dof, Sigma));
  EXPECT_THROW(wishart_log(Y_non_pos_def, dof, Sigma), std::domain_error);
  EXPECT_THROW(wishart_log(Y, dof, Sigma_non_pos_def), std::domain_error);
}

TEST(ProbDistributionsWishart, 2x2) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, Dynamic> Sigma(2, 2);
  Sigma << 1.848220, 1.899623, 1.899623, 12.751941;

  Matrix<double, Dynamic, Dynamic> Y(2, 2);
  Y << 2.011108, -11.206611, -11.206611, 112.94139;

  unsigned int dof = 3;

  // computed with MCMCpack in R
  double lp = log(8.658e-07);

  EXPECT_NEAR(lp, stan::math::wishart_log(Y, dof, Sigma), 0.01);
}
TEST(ProbDistributionsWishart, 4x4) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, Dynamic> Y(4, 4);
  Y << 7.988168, -9.555605, -14.47483, 4.395895, -9.555605, 44.750570, 49.21577,
      -18.454186, -14.474830, 49.21577, 60.08987, -21.48108, 4.395895,
      -18.454186, -21.48108, 7.885833;

  Matrix<double, Dynamic, Dynamic> Sigma(4, 4);
  Sigma << 2.9983662, 0.2898776, -2.650523, 0.1055911, 0.2898776, 11.4803610,
      7.1579931, -3.1129955, -2.650523, 7.1579931, 11.676181, -3.5866852,
      0.1055911, -3.1129955, -3.5866852, 1.4482736;

  double dof = 4;
  double log_p = log(8.034197e-10);
  EXPECT_NEAR(log_p, stan::math::wishart_log(Y, dof, Sigma), 0.01);

  dof = 5;
  log_p = log(1.517951e-10);
  EXPECT_NEAR(log_p, stan::math::wishart_log(Y, dof, Sigma), 0.01);
}

TEST(ProbDistributionsWishart, 2x2Propto) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, Dynamic> Sigma(2, 2);
  Sigma << 1.848220, 1.899623, 1.899623, 12.751941;

  Matrix<double, Dynamic, Dynamic> Y(2, 2);
  Y << 2.011108, -11.206611, -11.206611, 112.94139;

  unsigned int dof = 3;

  EXPECT_FLOAT_EQ(0.0, stan::math::wishart_log<true>(Y, dof, Sigma));
}

TEST(ProbDistributionsWishart, 4x4Propto) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, Dynamic> Y(4, 4);
  Y << 7.988168, -9.555605, -14.47483, 4.395895, -9.555605, 44.750570,
      49.215769, -18.454186, -14.474830, 49.215769, 60.08987, -21.48108,
      4.395895, -18.454186, -21.48108, 7.885833;

  Matrix<double, Dynamic, Dynamic> Sigma(4, 4);
  Sigma << 2.9983662, 0.2898776, -2.650523, 0.1055911, 0.2898776, 11.4803610,
      7.1579931, -3.1129955, -2.650523, 7.1579931, 11.676181, -3.5866852,
      0.1055911, -3.1129955, -3.5866852, 1.4482736;

  double dof = 4;
  EXPECT_FLOAT_EQ(0.0, stan::math::wishart_log<true>(Y, dof, Sigma));

  dof = 5;
  EXPECT_FLOAT_EQ(0.0, stan::math::wishart_log<true>(Y, dof, Sigma));
}

TEST(ProbDistributionsWishart, error) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::wishart_log;
  Matrix<double, Dynamic, Dynamic> Sigma;
  Matrix<double, Dynamic, Dynamic> Y;
  double nu;

  Sigma.resize(1, 1);
  Y.resize(1, 1);
  Sigma.setIdentity();
  Y.setIdentity();
  nu = 1;
  EXPECT_NO_THROW(wishart_log(Y, nu, Sigma));

  nu = 5;
  Sigma.resize(2, 1);
  EXPECT_THROW(wishart_log(Y, nu, Sigma), std::invalid_argument);

  nu = 5;
  Sigma.resize(2, 2);
  Y.resize(2, 1);
  EXPECT_THROW(wishart_log(Y, nu, Sigma), std::invalid_argument);

  nu = 5;
  Sigma.resize(2, 2);
  Y.resize(3, 3);
  EXPECT_THROW(wishart_log(Y, nu, Sigma), std::invalid_argument);

  Sigma.resize(3, 3);
  Sigma.setIdentity();
  Y.resize(3, 3);
  Y.setIdentity();
  nu = 3;
  EXPECT_NO_THROW(wishart_log(Y, nu, Sigma));
  nu = 2;
  EXPECT_THROW(wishart_log(Y, nu, Sigma), std::domain_error);
}
