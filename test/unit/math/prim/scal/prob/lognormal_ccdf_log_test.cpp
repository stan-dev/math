#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbLognormal, ccdf_log_matches_lccdf) {
  double y = 0.8;
  double mu = 1.1;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::lognormal_lccdf(y, mu, sigma)),
                  (stan::math::lognormal_ccdf_log(y, mu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::lognormal_lccdf<double, double, double>(y, mu, sigma)),
      (stan::math::lognormal_ccdf_log<double, double, double>(y, mu, sigma)));
}
