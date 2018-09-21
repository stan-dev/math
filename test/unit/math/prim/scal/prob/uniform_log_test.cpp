#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbUniform, log_matches_lpdf) {
  double y = 0.8;
  double alpha = 0.4;
  double beta = 2.3;

  EXPECT_FLOAT_EQ((stan::math::uniform_lpdf(y, alpha, beta)),
                  (stan::math::uniform_log(y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::uniform_lpdf<true>(y, alpha, beta)),
                  (stan::math::uniform_log<true>(y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::uniform_lpdf<false>(y, alpha, beta)),
                  (stan::math::uniform_log<false>(y, alpha, beta)));
  EXPECT_FLOAT_EQ(
      (stan::math::uniform_lpdf<true, double, double, double>(y, alpha, beta)),
      (stan::math::uniform_log<true, double, double, double>(y, alpha, beta)));
  EXPECT_FLOAT_EQ(
      (stan::math::uniform_lpdf<false, double, double, double>(y, alpha, beta)),
      (stan::math::uniform_log<false, double, double, double>(y, alpha, beta)));
  EXPECT_FLOAT_EQ(
      (stan::math::uniform_lpdf<double, double, double>(y, alpha, beta)),
      (stan::math::uniform_log<double, double, double>(y, alpha, beta)));
}
