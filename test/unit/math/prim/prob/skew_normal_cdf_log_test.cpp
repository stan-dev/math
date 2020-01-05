#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
TEST(ProbSkewNormal, cdf_log_matches_lcdf) {
  double y = 0.8;
  double mu = 2;
  double sigma = 2.3;
  double alpha = -3;

  EXPECT_FLOAT_EQ((stan::math::skew_normal_lcdf(y, mu, sigma, alpha)),
                  (stan::math::skew_normal_cdf_log(y, mu, sigma, alpha)));
  EXPECT_FLOAT_EQ(
      (stan::math::skew_normal_lcdf<double, double, double, double>(
          y, mu, sigma, alpha)),
      (stan::math::skew_normal_cdf_log<double, double, double, double>(
          y, mu, sigma, alpha)));
}
