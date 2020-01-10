#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbCauchy, cdf_log_matches_lcdf) {
  double y = 2;
  double mu = 1.5;
  double sigma = 0.5;

  EXPECT_FLOAT_EQ((stan::math::cauchy_lcdf(y, mu, sigma)),
                  (stan::math::cauchy_cdf_log(y, mu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::cauchy_lcdf<double, double, double>(y, mu, sigma)),
      (stan::math::cauchy_cdf_log<double, double, double>(y, mu, sigma)));
}
