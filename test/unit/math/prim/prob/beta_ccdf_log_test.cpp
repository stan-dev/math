#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbBeta, ccdf_log_matches_lccdf) {
  double y = 0.8;
  double alpha = 1.1;
  double beta = 2.3;

  EXPECT_FLOAT_EQ((stan::math::beta_lccdf(y, alpha, beta)),
                  (stan::math::beta_ccdf_log(y, alpha, beta)));
  EXPECT_FLOAT_EQ(
      (stan::math::beta_lccdf<double, double, double>(y, alpha, beta)),
      (stan::math::beta_ccdf_log<double, double, double>(y, alpha, beta)));
}
