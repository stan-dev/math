#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(ProbLkjCov, log_matches_lpmf) {
  unsigned int K = 4;
  Eigen::MatrixXd y(K, K);
  Eigen::VectorXd mu(K);
  Eigen::VectorXd sigma(K);
  y.setZero();
  y.diagonal().setOnes();
  mu << 0.1, 0.2, -1.2, 0.3;
  sigma << 1, 0.2, 4, 0.5;
  double eta = 1.2;

  EXPECT_FLOAT_EQ((stan::math::lkj_cov_lpdf(y, mu, sigma, eta)),
                  (stan::math::lkj_cov_log(y, mu, sigma, eta)));
  EXPECT_FLOAT_EQ((stan::math::lkj_cov_lpdf<true>(y, mu, sigma, eta)),
                  (stan::math::lkj_cov_log<true>(y, mu, sigma, eta)));
  EXPECT_FLOAT_EQ((stan::math::lkj_cov_lpdf<false>(y, mu, sigma, eta)),
                  (stan::math::lkj_cov_log<false>(y, mu, sigma, eta)));
  EXPECT_FLOAT_EQ(
      (stan::math::lkj_cov_lpdf<true, double, double>(y, mu, sigma, eta)),
      (stan::math::lkj_cov_log<true, double, double>(y, mu, sigma, eta)));
  EXPECT_FLOAT_EQ(
      (stan::math::lkj_cov_lpdf<false, double, double>(y, mu, sigma, eta)),
      (stan::math::lkj_cov_log<false, double, double>(y, mu, sigma, eta)));
}
