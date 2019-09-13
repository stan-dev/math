#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbBinomialLogit, log_matches_lpmf) {
  int n = 3;
  int N = 6;
  double theta = 1.2;

  EXPECT_FLOAT_EQ((stan::math::binomial_logit_lpmf(n, N, theta)),
                  (stan::math::binomial_logit_log(n, N, theta)));
  EXPECT_FLOAT_EQ((stan::math::binomial_logit_lpmf<true>(n, N, theta)),
                  (stan::math::binomial_logit_log<true>(n, N, theta)));
  EXPECT_FLOAT_EQ((stan::math::binomial_logit_lpmf<false>(n, N, theta)),
                  (stan::math::binomial_logit_log<false>(n, N, theta)));
  EXPECT_FLOAT_EQ((stan::math::binomial_logit_lpmf<true>(n, N, theta)),
                  (stan::math::binomial_logit_log<true>(n, N, theta)));
  EXPECT_FLOAT_EQ((stan::math::binomial_logit_lpmf<false>(n, N, theta)),
                  (stan::math::binomial_logit_log<false>(n, N, theta)));
  EXPECT_FLOAT_EQ((stan::math::binomial_logit_lpmf(n, N, theta)),
                  (stan::math::binomial_logit_log(n, N, theta)));
}
