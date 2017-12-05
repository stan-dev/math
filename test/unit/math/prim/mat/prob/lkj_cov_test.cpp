#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <random>

TEST(ProbDistributionsLkjCorr, testIdentity) {
  std::random_device rd;
  std::mt19937 mt(rd());
  unsigned int K = 4;
  Eigen::MatrixXd Sigma(K, K);
  Sigma.setZero();
  Sigma.diagonal().setOnes();
  double eta = rd() / static_cast<double>(RAND_MAX) + 0.5;
  double f = stan::math::do_lkj_constant(eta, K);
  EXPECT_FLOAT_EQ(f, stan::math::lkj_corr_log(Sigma, eta));
  eta = 1.0;
  f = stan::math::do_lkj_constant(eta, K);
  EXPECT_FLOAT_EQ(f, stan::math::lkj_corr_log(Sigma, eta));
}

TEST(ProbDistributionsLkjCorr, testHalf) {
  std::random_device rd;
  std::mt19937 mt(rd());
  unsigned int K = 4;
  Eigen::MatrixXd Sigma(K, K);
  Sigma.setConstant(0.5);
  Sigma.diagonal().setOnes();
  double eta = rd() / static_cast<double>(RAND_MAX) + 0.5;
  double f = stan::math::do_lkj_constant(eta, K);
  EXPECT_FLOAT_EQ(f + (eta - 1.0) * log(0.3125),
                  stan::math::lkj_corr_log(Sigma, eta));
  eta = 1.0;
  f = stan::math::do_lkj_constant(eta, K);
  EXPECT_FLOAT_EQ(f, stan::math::lkj_corr_log(Sigma, eta));
}

TEST(ProbDistributionsLkjCorr, Sigma) {
  std::random_device rd;
  std::mt19937 mt(rd());
  unsigned int K = 4;
  Eigen::MatrixXd Sigma(K, K);
  Sigma.setZero();
  Sigma.diagonal().setOnes();
  double eta = rd() / static_cast<double>(RAND_MAX) + 0.5;
  EXPECT_NO_THROW(stan::math::lkj_corr_log(Sigma, eta));

  EXPECT_THROW(stan::math::lkj_corr_log(Sigma, -eta), std::domain_error);

  Sigma = Sigma * -1.0;
  EXPECT_THROW(stan::math::lkj_corr_log(Sigma, eta), std::domain_error);
  Sigma = Sigma * (0.0 / 0.0);
  EXPECT_THROW(stan::math::lkj_corr_log(Sigma, eta), std::domain_error);
}
