#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbWeibull, log_matches_lpdf) {
  double y = 0.8;
  double alpha = 1.1;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::weibull_lpdf(y, alpha, sigma)),
                  (stan::math::weibull_log(y, alpha, sigma)));
  EXPECT_FLOAT_EQ((stan::math::weibull_lpdf<true>(y, alpha, sigma)),
                  (stan::math::weibull_log<true>(y, alpha, sigma)));
  EXPECT_FLOAT_EQ((stan::math::weibull_lpdf<false>(y, alpha, sigma)),
                  (stan::math::weibull_log<false>(y, alpha, sigma)));
  EXPECT_FLOAT_EQ((stan::math::weibull_lpdf<true>(y, alpha, sigma)),
                  (stan::math::weibull_log<true>(y, alpha, sigma)));
  EXPECT_FLOAT_EQ((stan::math::weibull_lpdf<false>(y, alpha, sigma)),
                  (stan::math::weibull_log<false>(y, alpha, sigma)));
  EXPECT_FLOAT_EQ((stan::math::weibull_lpdf(y, alpha, sigma)),
                  (stan::math::weibull_log(y, alpha, sigma)));
}
