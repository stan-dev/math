#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(ProbLkjCorrCholesky, log_matches_lpmf) {
  unsigned int K = 4;
  Eigen::MatrixXd Sigma(K,K);
  double eta = 1.2;
  Sigma.setZero();
  Sigma.diagonal().setOnes();
  EXPECT_FLOAT_EQ((stan::math::lkj_corr_cholesky_lpdf(Sigma, eta)),
                  (stan::math::lkj_corr_cholesky_log(Sigma, eta)));
  EXPECT_FLOAT_EQ((stan::math::lkj_corr_cholesky_lpdf<true>(Sigma, eta)),
                  (stan::math::lkj_corr_cholesky_log<true>(Sigma, eta)));
  EXPECT_FLOAT_EQ((stan::math::lkj_corr_cholesky_lpdf<false>(Sigma, eta)),
                  (stan::math::lkj_corr_cholesky_log<false>(Sigma, eta)));
  EXPECT_FLOAT_EQ((stan::math::lkj_corr_cholesky_lpdf<true, double, double>(Sigma, eta)),
                  (stan::math::lkj_corr_cholesky_log<true, double, double>(Sigma, eta)));
  EXPECT_FLOAT_EQ((stan::math::lkj_corr_cholesky_lpdf<false, double, double>(Sigma, eta)),
                  (stan::math::lkj_corr_cholesky_log<false, double, double>(Sigma, eta)));
}
