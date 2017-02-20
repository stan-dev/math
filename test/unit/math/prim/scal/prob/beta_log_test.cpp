#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbBeta, log_matches_lpdf) {
  double y = 0.8;
  double alpha = 1.1;
  double beta = 2.3;

  EXPECT_FLOAT_EQ((stan::math::beta_lpdf(y, alpha, beta)),
                  (stan::math::beta_log(y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::beta_lpdf<true>(y, alpha, beta)),
                  (stan::math::beta_log<true>(y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::beta_lpdf<false>(y, alpha, beta)),
                  (stan::math::beta_log<false>(y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::beta_lpdf<true, double, double, double>(y, alpha, beta)),
                  (stan::math::beta_log<true, double, double, double>(y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::beta_lpdf<false, double, double, double>(y, alpha, beta)),
                  (stan::math::beta_log<false, double, double, double>(y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::beta_lpdf<double, double, double>(y, alpha, beta)),
                  (stan::math::beta_log<double, double, double>(y, alpha, beta)));
}
