#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbBernoulli, log_matches_lpmf) {
  int n = 1;
  double theta = 0.3;

  EXPECT_FLOAT_EQ((stan::math::bernoulli_lpmf(n, theta)),
                  (stan::math::bernoulli_log(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_lpmf<true>(n, theta)),
                  (stan::math::bernoulli_log<true>(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_lpmf<false>(n, theta)),
                  (stan::math::bernoulli_log<false>(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_lpmf<true, double>(n, theta)),
                  (stan::math::bernoulli_log<true, double>(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_lpmf<false, double>(n, theta)),
                  (stan::math::bernoulli_log<false, double>(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_lpmf<double>(n, theta)),
                  (stan::math::bernoulli_log<double>(n, theta)));
}
