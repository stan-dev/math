#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbWeibull, ccdf_log_matches_lccdf) {
  double y = 0.8;
  double alpha = 1.1;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::weibull_lccdf(y, alpha, sigma)),
                  (stan::math::weibull_ccdf_log(y, alpha, sigma)));
  EXPECT_FLOAT_EQ((stan::math::weibull_lccdf<double, double, double>(y, alpha, sigma)),
                  (stan::math::weibull_ccdf_log<double, double, double>(y, alpha, sigma)));
}
