#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <test/unit/math/prim/scal/prob/util.hpp>
#include <vector>

TEST(ProbDistributionsUniform, error_check) {
  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::uniform_rng(1.0, 2.0, rng));

  EXPECT_THROW(
      stan::math::uniform_rng(stan::math::negative_infinity(), 2.0, rng),
      std::domain_error);
  EXPECT_THROW(stan::math::uniform_rng(1, stan::math::positive_infinity(), rng),
               std::domain_error);
  EXPECT_THROW(stan::math::uniform_rng(2, 1, rng), std::domain_error);
}

TEST(ProbDistributionsUniform, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000;
  int K = stan::math::round(2 * std::pow(N, 0.4));

  // Generate some samples.
  std::vector<double> samples;
  for (int i = 0; i < N; ++i) {
    samples.push_back(stan::math::uniform_rng(1.0, 2.0, rng));
  }

  // Generate quantiles for the uniform distribution.
  std::vector<double> quantiles;
  for (int i = 1; i <= K; ++i) {
    double frac = static_cast<double>(i) / K;
    quantiles.push_back(1.0 + frac);
  }

  // Assert that they match.
  assert_matches_quantiles(samples, quantiles, 1e-6);
}
