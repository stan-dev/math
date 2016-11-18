#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbWeibull, cdf_log_matches_lcdf) {
  double y = 0.8;
  double alpha = 1.1;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::weibull_lcdf(y, alpha, sigma)),
                  (stan::math::weibull_cdf_log(y, alpha, sigma)));
  EXPECT_FLOAT_EQ((stan::math::weibull_lcdf<double, double, double>(y, alpha, sigma)),
                  (stan::math::weibull_cdf_log<double, double, double>(y, alpha, sigma)));
}
