#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbBernoulliLogit, log_matches_lpmf) {
  int n = 1;
  double theta = 1.2;

  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf(n, theta)),
                  (stan::math::bernoulli_logit_log(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<true>(n, theta)),
                  (stan::math::bernoulli_logit_log<true>(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<false>(n, theta)),
                  (stan::math::bernoulli_logit_log<false>(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<true>(n, theta)),
                  (stan::math::bernoulli_logit_log<true>(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<false>(n, theta)),
                  (stan::math::bernoulli_logit_log<false>(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf(n, theta)),
                  (stan::math::bernoulli_logit_log(n, theta)));
}
