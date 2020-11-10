#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>

TEST(ProbDistributionsLkjCov, testIdentity) {
  boost::random::mt19937 rng;
  unsigned int K = 4;
  Eigen::MatrixXd Sigma(K, K);
  Sigma.setZero();
  Sigma.diagonal().setOnes();
  double mu = 0;
  Eigen::VectorXd muV = Eigen::VectorXd::Zero(K);
  double sd = 1;
  Eigen::VectorXd sdV = Eigen::VectorXd::Ones(K);
  double eta = stan::math::uniform_rng(0.5, 1.5, rng);
  double f = stan::math::do_lkj_constant(eta, K)
             + K * stan::math::lognormal_lpdf(1, 0, 1);
  EXPECT_FLOAT_EQ(f, stan::math::lkj_cov_lpdf(Sigma, mu, sd, eta));
  EXPECT_FLOAT_EQ(f, stan::math::lkj_cov_lpdf(Sigma, muV, sdV, eta));
  eta = 1.0;
  f = stan::math::do_lkj_constant(eta, K)
      + K * stan::math::lognormal_lpdf(1, 0, 1);
  EXPECT_FLOAT_EQ(f, stan::math::lkj_cov_lpdf(Sigma, mu, sd, eta));
  EXPECT_FLOAT_EQ(f, stan::math::lkj_cov_lpdf(Sigma, muV, sdV, eta));
}

TEST(ProbDistributionsLkjCov, testHalf) {
  boost::random::mt19937 rng;
  unsigned int K = 4;
  Eigen::MatrixXd Sigma(K, K);
  Sigma.setConstant(0.5);
  Sigma.diagonal().setOnes();
  double mu = 0;
  Eigen::VectorXd muV = Eigen::VectorXd::Zero(K);
  double sd = 1;
  Eigen::VectorXd sdV = Eigen::VectorXd::Ones(K);
  double eta = stan::math::uniform_rng(0.5, 1.5, rng);
  double f = stan::math::do_lkj_constant(eta, K)
             + K * stan::math::lognormal_lpdf(1, 0, 1)
             + (eta - 1.0) * log(0.3125);
  EXPECT_FLOAT_EQ(f, stan::math::lkj_cov_lpdf(Sigma, mu, sd, eta));
  EXPECT_FLOAT_EQ(f, stan::math::lkj_cov_lpdf(Sigma, muV, sdV, eta));
  eta = 1.0;
  f = stan::math::do_lkj_constant(eta, K)
      + K * stan::math::lognormal_lpdf(1, 0, 1);
  EXPECT_FLOAT_EQ(f, stan::math::lkj_cov_lpdf(Sigma, mu, sd, eta));
  EXPECT_FLOAT_EQ(f, stan::math::lkj_cov_lpdf(Sigma, muV, sdV, eta));
}

TEST(ProbDistributionsLkjCov, ErrorChecks) {
  boost::random::mt19937 rng;
  unsigned int K = 4;
  Eigen::MatrixXd Sigma(K, K);
  Sigma.setZero();
  Sigma.diagonal().setOnes();
  double mu = 0;
  Eigen::VectorXd muV = Eigen::VectorXd::Zero(K);
  double sd = 1;
  Eigen::VectorXd sdV = Eigen::VectorXd::Ones(K);
  double eta = stan::math::uniform_rng(0.5, 1.5, rng);
  EXPECT_NO_THROW(stan::math::lkj_cov_lpdf(Sigma, mu, sd, eta));

  // Error checks for non vectorized version.

  EXPECT_NO_THROW(stan::math::lkj_cov_lpdf(Sigma, -1, sd, eta));
  EXPECT_THROW(
      stan::math::lkj_cov_lpdf(Sigma, stan::math::NOT_A_NUMBER, sd, eta),
      std::domain_error);
  EXPECT_THROW(stan::math::lkj_cov_lpdf(Sigma, stan::math::INFTY, sd, eta),
               std::domain_error);
  EXPECT_THROW(
      stan::math::lkj_cov_lpdf(Sigma, stan::math::NEGATIVE_INFTY, sd, eta),
      std::domain_error);

  EXPECT_THROW(stan::math::lkj_cov_lpdf(Sigma, mu, -1, eta), std::domain_error);
  EXPECT_THROW(
      stan::math::lkj_cov_lpdf(Sigma, mu, stan::math::NOT_A_NUMBER, eta),
      std::domain_error);
  EXPECT_THROW(stan::math::lkj_cov_lpdf(Sigma, mu, stan::math::INFTY, eta),
               std::domain_error);
  EXPECT_THROW(
      stan::math::lkj_cov_lpdf(Sigma, mu, stan::math::NEGATIVE_INFTY, eta),
      std::domain_error);

  EXPECT_THROW(stan::math::lkj_cov_lpdf(Sigma, mu, sd, -1), std::domain_error);
  EXPECT_THROW(
      stan::math::lkj_cov_lpdf(Sigma, mu, sd, stan::math::NOT_A_NUMBER),
      std::domain_error);
  EXPECT_NO_THROW(stan::math::lkj_cov_lpdf(Sigma, mu, sd, stan::math::INFTY));
  EXPECT_THROW(
      stan::math::lkj_cov_lpdf(Sigma, mu, sd, stan::math::NEGATIVE_INFTY),
      std::domain_error);

  // Vectorized.

  Eigen::VectorXd muV1 = -muV;
  EXPECT_NO_THROW(stan::math::lkj_cov_lpdf(Sigma, muV1, sdV, eta));
  muV1 = stan::math::NOT_A_NUMBER * muV;
  EXPECT_THROW(stan::math::lkj_cov_lpdf(Sigma, muV1, sdV, eta),
               std::domain_error);
  muV1 = stan::math::INFTY * muV;
  EXPECT_THROW(stan::math::lkj_cov_lpdf(Sigma, muV1, sdV, eta),
               std::domain_error);
  muV1 = stan::math::NEGATIVE_INFTY * muV;
  EXPECT_THROW(stan::math::lkj_cov_lpdf(Sigma, muV1, sdV, eta),
               std::domain_error);

  Eigen::VectorXd sdV1 = -sdV;
  EXPECT_THROW(stan::math::lkj_cov_lpdf(Sigma, muV, sdV1, eta),
               std::domain_error);
  sdV1 = stan::math::NOT_A_NUMBER * sdV;
  EXPECT_THROW(stan::math::lkj_cov_lpdf(Sigma, muV, sdV1, eta),
               std::domain_error);
  sdV1 = stan::math::INFTY * sdV;
  EXPECT_THROW(stan::math::lkj_cov_lpdf(Sigma, muV, sdV1, eta),
               std::domain_error);
  sdV1 = stan::math::NEGATIVE_INFTY * sdV;
  EXPECT_THROW(stan::math::lkj_cov_lpdf(Sigma, muV, sdV1, eta),
               std::domain_error);

  EXPECT_THROW(stan::math::lkj_cov_lpdf(Sigma, muV, sdV, -1),
               std::domain_error);
  EXPECT_THROW(
      stan::math::lkj_cov_lpdf(Sigma, muV, sdV, stan::math::NOT_A_NUMBER),
      std::domain_error);
  EXPECT_NO_THROW(stan::math::lkj_cov_lpdf(Sigma, muV, sdV, stan::math::INFTY));
  EXPECT_THROW(
      stan::math::lkj_cov_lpdf(Sigma, muV, sdV, stan::math::NEGATIVE_INFTY),
      std::domain_error);
}
