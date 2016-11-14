#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbBetaBinomial, cdf_log_matches_lcdf) {
  int n = 2;
  int N = 6;
  double alpha = 1.1;
  double beta = 0.3;

  EXPECT_FLOAT_EQ((stan::math::beta_binomial_lcdf(n, N, alpha, beta)),
                  (stan::math::beta_binomial_cdf_log(n, N, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::beta_binomial_lcdf<int, int, double, double>(n, N, alpha, beta)),
                  (stan::math::beta_binomial_cdf_log<int, int, double, double>(n, N, alpha, beta)));
}
