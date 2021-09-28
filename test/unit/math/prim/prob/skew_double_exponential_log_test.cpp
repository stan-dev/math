#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbDistributionsSkewedDoubleExponential, log_matches_lpdf) {
  double y = 0.8;
  double mu = 2;
  double sigma = 2.3;
  double tau = .25;

  EXPECT_FLOAT_EQ((stan::math::skew_double_exponential_lpdf(y, mu, sigma, tau)),
                  (stan::math::skew_double_exponential_log(y, mu, sigma, tau)));
  EXPECT_FLOAT_EQ(
      (stan::math::skew_double_exponential_lpdf<true>(y, mu, sigma, tau)),
      (stan::math::skew_double_exponential_log<true>(y, mu, sigma, tau)));
  EXPECT_FLOAT_EQ(
      (stan::math::skew_double_exponential_lpdf<false>(y, mu, sigma, tau)),
      (stan::math::skew_double_exponential_log<false>(y, mu, sigma, tau)));
  EXPECT_FLOAT_EQ(
      (stan::math::skew_double_exponential_lpdf<true, double, double, double,
                                                double>(y, mu, sigma, tau)),
      (stan::math::skew_double_exponential_log<true, double, double, double,
                                               double>(y, mu, sigma, tau)));
  EXPECT_FLOAT_EQ(
      (stan::math::skew_double_exponential_lpdf<false, double, double, double,
                                                double>(y, mu, sigma, tau)),
      (stan::math::skew_double_exponential_log<false, double, double, double,
                                               double>(y, mu, sigma, tau)));
  EXPECT_FLOAT_EQ(
      (stan::math::skew_double_exponential_lpdf<double, double, double, double>(
          y, mu, sigma, tau)),
      (stan::math::skew_double_exponential_log<double, double, double, double>(
          y, mu, sigma, tau)));
}
