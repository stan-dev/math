#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbLogistic, ccdf_log_matches_lccdf) {
  double y = 0.8;
  double mu = -1.2;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::logistic_lccdf(y, mu, sigma)),
                  (stan::math::logistic_ccdf_log(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::logistic_lccdf<double, double, double>(y, mu, sigma)),
                  (stan::math::logistic_ccdf_log<double, double, double>(y, mu, sigma)));
}
