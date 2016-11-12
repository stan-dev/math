#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbBetaBinomial, log_matches_lpmf) {
  double n = 2;
  double N = 5.5;
  double alpha = 1.1;
  double beta = 0.3;

  EXPECT_FLOAT_EQ((stan::math::beta_binomial_lpmf(n, N, alpha, beta)),
                  (stan::math::beta_binomial_log(n, N, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::beta_binomial_lpmf<true>(n, N, alpha, beta)),
                  (stan::math::beta_binomial_log<true>(n, N, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::beta_binomial_lpmf<false>(n, N, alpha, beta)),
                  (stan::math::beta_binomial_log<false>(n, N, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::beta_binomial_lpmf<true, double, double, double, double>(n, N, alpha, beta)),
                  (stan::math::beta_binomial_log<true, double, double, double, double>(n, N, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::beta_binomial_lpmf<false, double, double, double, double>(n, N, alpha, beta)),
                  (stan::math::beta_binomial_log<false, double, double, double, double>(n, N, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::beta_binomial_lpmf<double, double, double, double>(n, N, alpha, beta)),
                  (stan::math::beta_binomial_log<double, double, double, double>(n, N, alpha, beta)));
}
