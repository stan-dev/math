#include <stan/math/prim/mat.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/meta/max_size.hpp>
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

template<typename T1, typename T2>
void chiSquareTestNormal(const T1& mu,
                         const T2& sigma) {
  boost::random::mt19937 rng;
  int N = 10000;
  int M = stan::max_size(mu, sigma);
  std::vector<std::vector<double> > samples_to_test_transpose;

  for (int n = 0; n < N; ++n) {
    samples_to_test_transpose.push_back(stan::math::normal_rng(mu, sigma, rng));
  }

  stan::scalar_seq_view<T1> mu_vec(mu);
  stan::scalar_seq_view<T2> sigma_vec(sigma);
  for (int m = 0; m < M; ++m) {
    int K = boost::math::round(2 * std::pow(N, 0.4));

    boost::math::normal_distribution<> dist (mu_vec[m], sigma_vec[m]);
    std::vector<double> quantiles;
    for (int i = 1; i < K; ++i) {
      double frac = static_cast<double>(i) / K;
      quantiles.push_back(quantile(dist, frac));
    }
    quantiles.push_back(std::numeric_limits<double>::max());

    std::vector<double> samples_to_test;
    for (int n = 0; n < N; ++n) {
      samples_to_test.push_back(samples_to_test_transpose[n][m]);
    }
    assert_matches_quantiles(samples_to_test, quantiles, 1e-6);
  }
}

TEST(ProbDistributionsNormal, chiSquareGoodnessFitTest_vec) {
  int M = 10;

  std::vector<double> mu(M), sigma(M);
  Matrix<double, Dynamic, 1> muV(M), sigmaV(M);
  Matrix<double, 1, Dynamic> muRV(M), sigmaRV(M);
  for(int m = 0; m < M; ++m) {
    mu[m] = m + 1.0;
    sigma[m] = m + 2.0;
    muV[m] = mu[m];
    sigmaV[m] = sigma[m];
    muRV[m] = mu[m];
    sigmaRV[m] = sigma[m];
  }

  chiSquareTestNormal(mu, sigma);
  chiSquareTestNormal(1.5, sigma);
  chiSquareTestNormal(mu, 1.0);

  chiSquareTestNormal(muV, sigmaV);
  chiSquareTestNormal(1.5, sigmaV);
  chiSquareTestNormal(muV, 1.0);

  chiSquareTestNormal(muRV, sigmaRV);
  chiSquareTestNormal(1.5, sigmaRV);
  chiSquareTestNormal(muRV, 1.0);

  chiSquareTestNormal(mu, sigmaV);
  chiSquareTestNormal(muV, sigma);

  chiSquareTestNormal(mu, sigmaRV);
  chiSquareTestNormal(muRV, sigma);
  
  chiSquareTestNormal(muV, sigmaRV);
  chiSquareTestNormal(muRV, sigmaV);
}
