#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <cmath>

TEST(ProbBetaBinomial, cdf_log_matches_lcdf) {
  int n = 2;
  int N = 6;
  double alpha = 1.1;
  double beta = 0.3;

  EXPECT_NEAR(stan::math::beta_binomial_lcdf(n, N, alpha, beta),
    stan::math::beta_binomial_cdf_log(n, N, alpha, beta), 1e-8);
  EXPECT_NEAR((stan::math::beta_binomial_lcdf(n, N, alpha, beta)),
    stan::math::beta_binomial_cdf_log(n, N, alpha, beta), 1e-8);
}


TEST(ProbBetaBinomial, lcdf_like_lcdf) {
  int n = 10;
  int N = 10;
  double alpha = 3.0;
  double beta = 1.0;

  EXPECT_NEAR(0.0, stan::math::beta_binomial_lcdf(n, N, alpha, beta), 1e-8);
  EXPECT_NEAR(0.0, std::exp(stan::math::beta_binomial_lcdf(0.0, N, alpha, beta)), 1e-8);

}

TEST(ProbBetaBinomial, lcdf_matches_lpmf) {
  int n = 9;
  int N = 10;
  double alpha = 3.0;
  double beta = 2.1;

  double pmf_sum = 0.0;
  for (int i = 0; i <= n; ++i) 
    pmf_sum += std::exp(stan::math::beta_binomial_lpmf(i, N, alpha, beta));

  EXPECT_NEAR(
    pmf_sum,
    std::exp(stan::math::beta_binomial_lcdf(n, N, alpha, beta)),
    1e-8);
}

TEST(ProbBetaBinomial, lcdf_matches_mathematica) {
  int n = 8;
  int N = 10;
  double alpha = 3.0;
  double beta = 1.0;

  //  EXPECT_NEAR(-0.5500463, (stan::math::beta_binomial_lcdf(n, N, alpha, beta)), 1e-8);
  // FIXME: this point _should_ be defined for the beta_binomial_lcdf to be defined 
  // over its full parameter range but the power-series is not defined.  Leaving the test
  // in place with the current behavior.
  EXPECT_THROW(stan::math::beta_binomial_lcdf(n, N, alpha, beta), std::domain_error);
}


