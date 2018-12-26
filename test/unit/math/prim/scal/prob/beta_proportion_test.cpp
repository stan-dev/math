#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <test/unit/math/prim/scal/prob/util.hpp>
#include <limits>
#include <vector>

TEST(ProbDistributionsBetaProportion, error_check) {
  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::beta_proportion_rng(0.5, 3.0, rng));
  EXPECT_NO_THROW(stan::math::beta_proportion_rng(1e-10, 1e-10, rng));
  EXPECT_THROW(stan::math::beta_proportion_rng(2.0, -1.0, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::beta_proportion_rng(-2.0, 1.0, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::beta_proportion_rng(-2.0, -1.0, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::beta_proportion_rng(stan::math::positive_infinity(),
                                               1.0, rng),
               std::domain_error);
  EXPECT_THROW(
      stan::math::beta_proportion_rng(2, stan::math::positive_infinity(), rng),
      std::domain_error);
}

TEST(ProbDistributionsBetaProportion, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000;
  int K = stan::math::round(2 * std::pow(N, 0.4));

  double mu = 0.5;     // location
  double kappa = 3.0;  // precision

  std::vector<double> samples;
  for (int i = 0; i < N; ++i) {
    samples.push_back(stan::math::beta_proportion_rng(mu, kappa, rng));
  }

  // transform from location and precision parameterization
  // into shape1 (alpha) and shape2 (beta) parameterization
  double alpha = mu * kappa;
  double beta = kappa - alpha;

  // Generate quantiles from boost's beta distribution
  boost::math::beta_distribution<> dist(alpha, beta);
  std::vector<double> quantiles;
  for (int i = 1; i < K; ++i) {
    double frac = static_cast<double>(i) / K;
    quantiles.push_back(quantile(dist, frac));
  }
  quantiles.push_back(std::numeric_limits<double>::max());

  // Assert that they match
  assert_matches_quantiles(samples, quantiles, 1e-6);
}

TEST(ProbDistributionsBetaProportion, chiSquareGoodnessFitTest2) {
  boost::random::mt19937 rng;
  int N = 10000;
  int K = stan::math::round(2 * std::pow(N, 0.4));

  double mu = 0.3;     // location
  double kappa = 0.5;  // precision

  std::vector<double> samples;
  for (int i = 0; i < N; ++i) {
    samples.push_back(stan::math::beta_proportion_rng(mu, kappa, rng));
  }

  // transform from location and precision parameterization
  // into shape1 (alpha) and shape2 (beta) parameterization
  double alpha = mu * kappa;
  double beta = kappa - alpha;

  // Generate quantiles from boost's beta distribution
  boost::math::beta_distribution<> dist(alpha, beta);
  std::vector<double> quantiles;
  for (int i = 1; i < K; ++i) {
    double frac = static_cast<double>(i) / K;
    quantiles.push_back(quantile(dist, frac));
  }
  quantiles.push_back(std::numeric_limits<double>::max());

  // Assert that they match
  assert_matches_quantiles(samples, quantiles, 1e-6);
}
