#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbBetaProportion, cdf_log_matches_lcdf) {
  double y = 0.8;
  double mu = .51;
  double kappa = 2.3;

  EXPECT_FLOAT_EQ((stan::math::beta_proportion_lcdf(y, mu, kappa)),
                  (stan::math::beta_proportion_cdf_log(y, mu, kappa)));
  EXPECT_FLOAT_EQ(
      (stan::math::beta_proportion_lcdf<double, double, double>(y, mu, kappa)),
      (stan::math::beta_proportion_cdf_log<double, double, double>(y, mu,
                                                                   kappa)));
}
