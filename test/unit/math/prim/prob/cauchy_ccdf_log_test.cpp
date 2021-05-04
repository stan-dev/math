#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbCauchy, ccdf_log_matches_lcdf) {
  double y = 2;
  double mu = 1.5;
  double sigma = 0.5;

  EXPECT_FLOAT_EQ((stan::math::cauchy_lccdf(y, mu, sigma)),
                  (stan::math::cauchy_ccdf_log(y, mu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::cauchy_lccdf<double, double, double>(y, mu, sigma)),
      (stan::math::cauchy_ccdf_log<double, double, double>(y, mu, sigma)));
}
