#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbLognormal, cdf_log_matches_lcdf) {
  double y = 0.8;
  double mu = 1.1;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::lognormal_lcdf(y, mu, sigma)),
                  (stan::math::lognormal_cdf_log(y, mu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::lognormal_lcdf<double, double, double>(y, mu, sigma)),
      (stan::math::lognormal_cdf_log<double, double, double>(y, mu, sigma)));
}
