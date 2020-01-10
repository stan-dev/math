#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbSkewNormal, ccdf_log_matches_lccdf) {
  double y = 0.8;
  double mu = 2;
  double sigma = 2.3;
  double alpha = -3;

  EXPECT_FLOAT_EQ((stan::math::skew_normal_lccdf(y, mu, sigma, alpha)),
                  (stan::math::skew_normal_ccdf_log(y, mu, sigma, alpha)));
  EXPECT_FLOAT_EQ(
      (stan::math::skew_normal_lccdf<double, double, double, double>(
          y, mu, sigma, alpha)),
      (stan::math::skew_normal_ccdf_log<double, double, double, double>(
          y, mu, sigma, alpha)));
}
