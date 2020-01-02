#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbExponential, log_matches_lpdf) {
  double y = 0.8;
  double beta = 1.2;

  EXPECT_FLOAT_EQ((stan::math::exponential_lpdf(y, beta)),
                  (stan::math::exponential_log(y, beta)));
  EXPECT_FLOAT_EQ((stan::math::exponential_lpdf<true>(y, beta)),
                  (stan::math::exponential_log<true>(y, beta)));
  EXPECT_FLOAT_EQ((stan::math::exponential_lpdf<false>(y, beta)),
                  (stan::math::exponential_log<false>(y, beta)));
  EXPECT_FLOAT_EQ((stan::math::exponential_lpdf<true, double, double>(y, beta)),
                  (stan::math::exponential_log<true, double, double>(y, beta)));
  EXPECT_FLOAT_EQ(
      (stan::math::exponential_lpdf<false, double, double>(y, beta)),
      (stan::math::exponential_log<false, double, double>(y, beta)));
  EXPECT_FLOAT_EQ((stan::math::exponential_lpdf<double, double>(y, beta)),
                  (stan::math::exponential_log<double, double>(y, beta)));
}
