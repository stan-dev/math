#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <cmath>

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


TEST(ProbBetaBinomial, lcdf_like_lcdf) {
  int n = 10;
  int N = 10;
  double alpha = 3.0;
  double beta = 1.0;

  EXPECT_FLOAT_EQ(0.0, (stan::math::beta_binomial_lcdf(n, N, alpha, beta)));
  EXPECT_FLOAT_EQ(0.0, std::exp(stan::math::beta_binomial_lcdf(0.0, N, alpha, beta)));
}

TEST(ProbBetaBinomial, lcdf_matches_mathematica) {
  int n = 8;
  int N = 10;
  double alpha = 3.0;
  double beta = 1.0;

  EXPECT_FLOAT_EQ(-0.5500463, (stan::math::beta_binomial_lcdf(n, N, alpha, beta)));
}


