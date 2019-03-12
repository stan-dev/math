#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <test/unit/math/prim/scal/prob/util.hpp>
#include <limits>
#include <vector>

TEST(ProbDistributionsFrechet, error_check) {
  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::frechet_rng(2.0, 3.0, rng));

  EXPECT_THROW(stan::math::frechet_rng(-2.0, 3.0, rng), std::domain_error);
  EXPECT_THROW(stan::math::frechet_rng(2.0, -3.0, rng), std::domain_error);
  EXPECT_THROW(
      stan::math::frechet_rng(stan::math::positive_infinity(), 3.0, rng),
      std::domain_error);
}

TEST(ProbDistributionsFrechet, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000;
  int K = stan::math::round(2 * std::pow(N, 0.4));

  // and invert them to compare to Weibull
  std::vector<double> samples;
  for (int i = 0; i < N; ++i) {
    samples.push_back(1.0 / stan::math::frechet_rng(2.0, 1.0 / 5.0, rng));
  }

  // Generate quantiles from boost's Frechet distribution
  // Actually, boost doesn't have a Frechet distribution
  // So we use the Weibull.
  boost::math::weibull_distribution<> dist(2.0, 5.0);
  std::vector<double> quantiles;
  for (int i = 1; i < K; ++i) {
    double frac = static_cast<double>(i) / K;
    quantiles.push_back(quantile(dist, frac));
  }
  quantiles.push_back(std::numeric_limits<double>::max());

  // Assert that they match
  assert_matches_quantiles(samples, quantiles, 1e-6);
}
