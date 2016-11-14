#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbDoubleExponential, log_matches_lpdf) {
  double y = 0.8;
  double mu = 1.1;
  double sigma = 2.3;
  
  EXPECT_FLOAT_EQ((stan::math::double_exponential_lpdf(y, mu, sigma)),
                  (stan::math::double_exponential_log(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::double_exponential_lpdf<true>(y, mu, sigma)),
                  (stan::math::double_exponential_log<true>(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::double_exponential_lpdf<false>(y, mu, sigma)),
                  (stan::math::double_exponential_log<false>(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::double_exponential_lpdf<true, double, double, double>(y, mu, sigma)),
                  (stan::math::double_exponential_log<true, double, double, double>(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::double_exponential_lpdf<false, double, double, double>(y, mu, sigma)),
                  (stan::math::double_exponential_log<false, double, double, double>(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::double_exponential_lpdf<double, double, double>(y, mu, sigma)),
                  (stan::math::double_exponential_log<double, double, double>(y, mu, sigma)));
}
