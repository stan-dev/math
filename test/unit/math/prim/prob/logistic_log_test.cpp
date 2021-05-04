#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbLogistic, log_matches_lpdf) {
  double y = 0.8;
  double mu = -1.2;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::logistic_lpdf(y, mu, sigma)),
                  (stan::math::logistic_log(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::logistic_lpdf<true>(y, mu, sigma)),
                  (stan::math::logistic_log<true>(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::logistic_lpdf<false>(y, mu, sigma)),
                  (stan::math::logistic_log<false>(y, mu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::logistic_lpdf<true, double, double, double>(y, mu, sigma)),
      (stan::math::logistic_log<true, double, double, double>(y, mu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::logistic_lpdf<false, double, double, double>(y, mu, sigma)),
      (stan::math::logistic_log<false, double, double, double>(y, mu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::logistic_lpdf<double, double, double>(y, mu, sigma)),
      (stan::math::logistic_log<double, double, double>(y, mu, sigma)));
}
