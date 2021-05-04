#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbBinomial, cdf_log_matches_lcdf) {
  int n = 3;
  int N = 6;
  double theta = 0.4;

  EXPECT_FLOAT_EQ((stan::math::binomial_lcdf(n, N, theta)),
                  (stan::math::binomial_cdf_log(n, N, theta)));
  EXPECT_FLOAT_EQ((stan::math::binomial_lcdf<double>(n, N, theta)),
                  (stan::math::binomial_cdf_log<double>(n, N, theta)));
}
