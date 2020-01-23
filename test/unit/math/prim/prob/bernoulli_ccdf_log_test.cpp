#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbBernoulli, ccdf_log_matches_lccdf) {
  int n = 1;
  double theta = 0.3;

  EXPECT_FLOAT_EQ((stan::math::bernoulli_lccdf(n, theta)),
                  (stan::math::bernoulli_ccdf_log(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_lccdf<double>(n, theta)),
                  (stan::math::bernoulli_ccdf_log<double>(n, theta)));
}
