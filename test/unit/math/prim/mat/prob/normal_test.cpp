#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <test/unit/math/prim/scal/prob/util.hpp>

using Eigen::Dynamic;
using Eigen::Matrix;

TEST(ProbDistributionsNormal, error_check_vec) {
  boost::random::mt19937 rng;

  double mu = 1.0;
  double sigma = 2.0;
  Matrix<double, Dynamic, 1> mu_vec(3), sigma_vec(3),
    mu_vec2(4), sigma_vec2(4);

  mu_vec << -1.0, 2.0, 3.0;
  sigma_vec << 1.0, 2.0, 3.0;
  mu_vec2 << 1.0, -2.0, 3.0, 4.0;
  sigma_vec2 << 1.0, 2.0, 3.0, 4.0;

  EXPECT_NO_THROW(stan::math::normal_rng(mu, sigma_vec, rng));
  EXPECT_NO_THROW(stan::math::normal_rng(mu_vec, sigma, rng));
  EXPECT_NO_THROW(stan::math::normal_rng(mu_vec, sigma_vec, rng));

  mu = stan::math::positive_infinity();
  mu_vec[1] = stan::math::positive_infinity();

  EXPECT_THROW(stan::math::normal_rng(mu, sigma_vec, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::normal_rng(mu_vec, sigma, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::normal_rng(mu_vec, sigma_vec, rng),
               std::domain_error);

  mu = stan::math::negative_infinity();
  mu_vec[1] = 1.0;
  mu_vec[0] = stan::math::negative_infinity();
  EXPECT_THROW(stan::math::normal_rng(mu, sigma_vec, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::normal_rng(mu_vec, sigma, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::normal_rng(mu_vec, sigma_vec, rng),
               std::domain_error);

  mu = 0.0;
  mu_vec[0] = 0.0;
  sigma = 0.0;
  sigma_vec[0] = 0.0;
  EXPECT_THROW(stan::math::normal_rng(mu, sigma_vec, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::normal_rng(mu_vec, sigma, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::normal_rng(mu_vec, sigma_vec, rng),
               std::domain_error);

  sigma = -2.0;
  sigma_vec[0] = 1.0;
  sigma_vec[1] = -2.0;
  EXPECT_THROW(stan::math::normal_rng(mu, sigma_vec, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::normal_rng(mu_vec, sigma, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::normal_rng(mu_vec, sigma_vec, rng),
               std::domain_error);

  sigma = 1.0;
  sigma_vec[1] = 1.0;
  EXPECT_THROW(stan::math::normal_rng(mu_vec, sigma_vec2, rng),
               std::invalid_argument);
  EXPECT_THROW(stan::math::normal_rng(mu_vec2, sigma_vec, rng),
  std::invalid_argument);
}

TEST(ProbDistributionsNormal, chiSquareGoodnessFitTest_vec) {
  boost::random::mt19937 rng;
  int N = 10000;
  int K = boost::math::round(2 * std::pow(N, 0.4));

  std::vector<double> samples(N);

  Matrix<double, Dynamic, 1> mu(N), sigma(N), r;
  for (int i = 0; i < N; ++i) {
    mu(i) = 2.0;
    sigma(i) = 1.0;
  }

  r = stan::math::normal_rng(mu, sigma, rng);

  for (int i = 0; i < N; ++i) {
    samples[i] = r(i);
  }

  //Generate quantiles from boost's normal distribution
  boost::math::normal_distribution<> dist (2.0, 1.0);
  std::vector<double> quantiles;
  for (int i = 1; i < K; ++i) {
    double frac = static_cast<double>(i) / K;
    quantiles.push_back(quantile(dist, frac));
  }
  quantiles.push_back(std::numeric_limits<double>::max());

  //Assert that they match
  assert_matches_quantiles(samples, quantiles, 1e-6);
}
