#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbBetaProportion, ccdf_log_matches_lccdf) {
  double y = 0.8;
  double mu = 0.51;
  double kappa = 2.3;

  EXPECT_FLOAT_EQ((stan::math::beta_proportion_lccdf(y, mu, kappa)),
                  (stan::math::beta_proportion_ccdf_log(y, mu, kappa)));
  EXPECT_FLOAT_EQ(
      (stan::math::beta_proportion_lccdf<double, double, double>(y, mu, kappa)),
      (stan::math::beta_proportion_ccdf_log<double, double, double>(y, mu,
                                                                    kappa)));
}
