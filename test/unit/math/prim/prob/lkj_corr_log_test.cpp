#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbLkjCorr, log_matches_lpmf) {
  unsigned int K = 4;
  Eigen::MatrixXd Sigma(K, K);
  Sigma.setZero();
  Sigma.diagonal().setOnes();
  double eta = 1.2;

  EXPECT_FLOAT_EQ((stan::math::lkj_corr_lpdf(Sigma, eta)),
                  (stan::math::lkj_corr_log(Sigma, eta)));
  EXPECT_FLOAT_EQ((stan::math::lkj_corr_lpdf<true>(Sigma, eta)),
                  (stan::math::lkj_corr_log<true>(Sigma, eta)));
  EXPECT_FLOAT_EQ((stan::math::lkj_corr_lpdf<false>(Sigma, eta)),
                  (stan::math::lkj_corr_log<false>(Sigma, eta)));
  EXPECT_FLOAT_EQ((stan::math::lkj_corr_lpdf<true>(Sigma, eta)),
                  (stan::math::lkj_corr_log<true>(Sigma, eta)));
  EXPECT_FLOAT_EQ((stan::math::lkj_corr_lpdf<false>(Sigma, eta)),
                  (stan::math::lkj_corr_log<false>(Sigma, eta)));
}
