#include <stan/math/prim/scal.hpp>
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
  EXPECT_FLOAT_EQ((stan::math::bernoulli_lpmf<true>(n, theta)),
                  (stan::math::bernoulli_log<true>(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_lpmf<false>(n, theta)),
                  (stan::math::bernoulli_log<false>(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_lpmf(n, theta)),
                  (stan::math::bernoulli_log(n, theta)));
}
