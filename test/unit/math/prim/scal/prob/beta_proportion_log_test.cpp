#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbBetaProportion, log_matches_lpdf) {
  double y = 0.8;
  double mu = .51;
  double kappa = 2.3;

  EXPECT_FLOAT_EQ((stan::math::beta_proportion_lpdf(y, mu, kappa)),
                  (stan::math::beta_proportion_log(y, mu, kappa)));
  EXPECT_FLOAT_EQ((stan::math::beta_proportion_lpdf<true>(y, mu, kappa)),
                  (stan::math::beta_proportion_log<true>(y, mu, kappa)));
  EXPECT_FLOAT_EQ((stan::math::beta_proportion_lpdf<false>(y, mu, kappa)),
                  (stan::math::beta_proportion_log<false>(y, mu, kappa)));
  EXPECT_FLOAT_EQ((stan::math::beta_proportion_lpdf<true>(y, mu, kappa)),
                  (stan::math::beta_proportion_log<true>(y, mu, kappa)));
  EXPECT_FLOAT_EQ((stan::math::beta_proportion_lpdf<false>(y, mu, kappa)),
                  (stan::math::beta_proportion_log<false>(y, mu, kappa)));
  EXPECT_FLOAT_EQ((stan::math::beta_proportion_lpdf(y, mu, kappa)),
                  (stan::math::beta_proportion_log(y, mu, kappa)));
}
