#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbBeta, cdf_log_matches_lcdf) {
  double y = 0.8;
  double alpha = 1.1;
  double beta = 2.3;

  EXPECT_FLOAT_EQ((stan::math::beta_lcdf(y, alpha, beta)),
                  (stan::math::beta_cdf_log(y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::beta_lcdf<double, double, double>(y, alpha, beta)),
                  (stan::math::beta_cdf_log<double, double, double>(y, alpha, beta)));
}
