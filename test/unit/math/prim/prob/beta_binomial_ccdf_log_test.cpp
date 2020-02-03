#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>

TEST(ProbBetaBinomial, ccdf_log_matches_lccdf) {
  int n = 2;
  int N = 6;
  double alpha = 1.1;
  double beta = 0.3;

  EXPECT_FLOAT_EQ((stan::math::beta_binomial_lccdf(n, N, alpha, beta)),
                  (stan::math::beta_binomial_ccdf_log(n, N, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::beta_binomial_lccdf<int, int, double, double>(
                      n, N, alpha, beta)),
                  (stan::math::beta_binomial_ccdf_log<int, int, double, double>(
                      n, N, alpha, beta)));
}

TEST(ProbBetaBinomial, lccdf_matches_lpmf) {
  int n = 10;
  int N = 10;
  double alpha = 4.0;
  double beta = 5.5;

  double pmf_sum = 1.0;
  for (int i = 0; i <= n; ++i) {
    pmf_sum -= std::exp(stan::math::beta_binomial_lpmf(i, N, alpha, beta));
    EXPECT_NEAR(pmf_sum,
                std::exp(stan::math::beta_binomial_lccdf(i, N, alpha, beta)),
                1e-8);
  }
}
