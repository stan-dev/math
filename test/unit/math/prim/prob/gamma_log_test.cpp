#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbGamma, log_matches_lpdf) {
  double y = 0.8;
  double alpha = 1.1;
  double beta = 2.3;

  EXPECT_FLOAT_EQ((stan::math::gamma_lpdf(y, alpha, beta)),
                  (stan::math::gamma_log(y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::gamma_lpdf<true>(y, alpha, beta)),
                  (stan::math::gamma_log<true>(y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::gamma_lpdf<false>(y, alpha, beta)),
                  (stan::math::gamma_log<false>(y, alpha, beta)));
  EXPECT_FLOAT_EQ(
      (stan::math::gamma_lpdf<true, double, double, double>(y, alpha, beta)),
      (stan::math::gamma_log<true, double, double, double>(y, alpha, beta)));
  EXPECT_FLOAT_EQ(
      (stan::math::gamma_lpdf<false, double, double, double>(y, alpha, beta)),
      (stan::math::gamma_log<false, double, double, double>(y, alpha, beta)));
  EXPECT_FLOAT_EQ(
      (stan::math::gamma_lpdf<double, double, double>(y, alpha, beta)),
      (stan::math::gamma_log<double, double, double>(y, alpha, beta)));
}
