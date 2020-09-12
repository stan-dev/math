#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbLognormal, log_matches_lpdf) {
  double y = 0.8;
  double mu = 1.1;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::lognormal_lpdf(y, mu, sigma)),
                  (stan::math::lognormal_log(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::lognormal_lpdf<true>(y, mu, sigma)),
                  (stan::math::lognormal_log<true>(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::lognormal_lpdf<false>(y, mu, sigma)),
                  (stan::math::lognormal_log<false>(y, mu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::lognormal_lpdf<true, double, double, double>(y, mu, sigma)),
      (stan::math::lognormal_log<true, double, double, double>(y, mu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::lognormal_lpdf<false, double, double, double>(y, mu, sigma)),
      (stan::math::lognormal_log<false, double, double, double>(y, mu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::lognormal_lpdf<double, double, double>(y, mu, sigma)),
      (stan::math::lognormal_log<double, double, double>(y, mu, sigma)));

  EXPECT_FLOAT_EQ((stan::math::lognormal_lpdf(0.8, 2.0, 2.0)),
                  (stan::math::lognormal_log(0.8, 2, 2)));
}
