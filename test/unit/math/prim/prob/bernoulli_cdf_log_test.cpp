#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbBernoulli, cdf_log_matches_lcdf) {
  int n = 1;
  double theta = 0.3;

  EXPECT_FLOAT_EQ((stan::math::bernoulli_lcdf(n, theta)),
                  (stan::math::bernoulli_cdf_log(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_lcdf<double>(n, theta)),
                  (stan::math::bernoulli_cdf_log<double>(n, theta)));
}
