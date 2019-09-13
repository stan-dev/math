#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbBinomial, log_matches_lpmf) {
  int n = 3;
  int N = 6;
  double theta = 0.4;

  EXPECT_FLOAT_EQ((stan::math::binomial_lpmf(n, N, theta)),
                  (stan::math::binomial_log(n, N, theta)));
  EXPECT_FLOAT_EQ((stan::math::binomial_lpmf<true>(n, N, theta)),
                  (stan::math::binomial_log<true>(n, N, theta)));
  EXPECT_FLOAT_EQ((stan::math::binomial_lpmf<false>(n, N, theta)),
                  (stan::math::binomial_log<false>(n, N, theta)));
  EXPECT_FLOAT_EQ((stan::math::binomial_lpmf<true>(n, N, theta)),
                  (stan::math::binomial_log<true>(n, N, theta)));
  EXPECT_FLOAT_EQ((stan::math::binomial_lpmf<false>(n, N, theta)),
                  (stan::math::binomial_log<false>(n, N, theta)));
  EXPECT_FLOAT_EQ((stan::math::binomial_lpmf(n, N, theta)),
                  (stan::math::binomial_log(n, N, theta)));
}
