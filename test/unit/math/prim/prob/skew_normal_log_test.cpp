#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbSkewNormal, log_matches_lpdf) {
  double y = 0.8;
  double mu = 2;
  double sigma = 2.3;
  double alpha = -3;

  EXPECT_FLOAT_EQ((stan::math::skew_normal_lpdf(y, mu, sigma, alpha)),
                  (stan::math::skew_normal_log(y, mu, sigma, alpha)));
  EXPECT_FLOAT_EQ((stan::math::skew_normal_lpdf<true>(y, mu, sigma, alpha)),
                  (stan::math::skew_normal_log<true>(y, mu, sigma, alpha)));
  EXPECT_FLOAT_EQ((stan::math::skew_normal_lpdf<false>(y, mu, sigma, alpha)),
                  (stan::math::skew_normal_log<false>(y, mu, sigma, alpha)));
  EXPECT_FLOAT_EQ(
      (stan::math::skew_normal_lpdf<true, double, double, double, double>(
          y, mu, sigma, alpha)),
      (stan::math::skew_normal_log<true, double, double, double, double>(
          y, mu, sigma, alpha)));
  EXPECT_FLOAT_EQ(
      (stan::math::skew_normal_lpdf<false, double, double, double, double>(
          y, mu, sigma, alpha)),
      (stan::math::skew_normal_log<false, double, double, double, double>(
          y, mu, sigma, alpha)));
  EXPECT_FLOAT_EQ((stan::math::skew_normal_lpdf<double, double, double, double>(
                      y, mu, sigma, alpha)),
                  (stan::math::skew_normal_log<double, double, double, double>(
                      y, mu, sigma, alpha)));
}
