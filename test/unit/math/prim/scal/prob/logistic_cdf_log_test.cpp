#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbLogistic, cdf_log_matches_lcdf) {
  double y = 0.8;
  double mu = -1.2;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::logistic_lcdf(y, mu, sigma)),
                  (stan::math::logistic_cdf_log(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::logistic_lcdf<double, double, double>(y, mu, sigma)),
                  (stan::math::logistic_cdf_log<double, double, double>(y, mu, sigma)));
}
