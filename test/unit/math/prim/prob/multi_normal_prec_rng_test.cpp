#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim.hpp>
#include <test/unit/math/prim/prob/util.hpp>
#include <limits>
#include <vector>

TEST(ProbDistributionsMultiNormalPrec, vectorized) {
  // Test scalar/vector combinations.
  boost::random::mt19937 rng;

  Eigen::VectorXd mu(3);
  Eigen::RowVectorXd mu_t(3);
  std::vector<Eigen::VectorXd> vec_mu(2);
  std::vector<Eigen::RowVectorXd> vec_mu_t(2);
  vec_mu[0] = mu;
  vec_mu_t[0] = mu;
  mu << 2.0, -1.0, 4.0;
  vec_mu[1] = mu;
  vec_mu_t[1] = mu;
  mu << 1.0, -1.0, 3.0;
  mu_t = mu;

  Eigen::MatrixXd Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;

  EXPECT_NO_THROW(stan::math::multi_normal_prec_rng(mu, Sigma, rng));
  EXPECT_NO_THROW(stan::math::multi_normal_prec_rng(mu_t, Sigma, rng));
  EXPECT_NO_THROW(stan::math::multi_normal_prec_rng(vec_mu, Sigma, rng));
  EXPECT_NO_THROW(stan::math::multi_normal_prec_rng(vec_mu_t, Sigma, rng));
}

TEST(ProbDistributionsMultiNormalPrec, policiesSigma) {
  boost::random::mt19937 rng;

  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();
  double ninf = -std::numeric_limits<double>::infinity();

  Eigen::VectorXd mu(2);
  mu << 1.0, -1.0;
  Eigen::MatrixXd Sigma(2, 2);
  Sigma << 9.0, -3.0, -3.0, 4.0;
  EXPECT_NO_THROW(stan::math::multi_normal_prec_rng(mu, Sigma, rng));

  // non-symmetric
  Sigma(0, 1) = -2.5;
  EXPECT_THROW(stan::math::multi_normal_prec_rng(mu, Sigma, rng),
               std::domain_error);
  Sigma(0, 1) = -3;

  // not positive-definite
  Sigma(0, 1) = -10;
  EXPECT_THROW(stan::math::multi_normal_prec_rng(mu, Sigma, rng),
               std::domain_error);
  Sigma(0, 1) = -3;

  // not square
  Eigen::MatrixXd rect_Sigma(2, 3);
  rect_Sigma << 1, 0, 0, 0, 1, 0;
  EXPECT_THROW(stan::math::multi_normal_prec_rng(mu, rect_Sigma, rng),
               std::invalid_argument);

  // NaN
  Sigma(0, 0) = nan;
  EXPECT_THROW(stan::math::multi_normal_prec_rng(mu, Sigma, rng),
               std::domain_error);
  Sigma(0, 0) = 9.0;

  // inf
  Sigma(0, 0) = inf;
  EXPECT_THROW(stan::math::multi_normal_prec_rng(mu, Sigma, rng),
               std::domain_error);
  Sigma(0, 0) = ninf;
  EXPECT_THROW(stan::math::multi_normal_prec_rng(mu, Sigma, rng),
               std::domain_error);
  Sigma(0, 0) = 9.0;
}

TEST(ProbDistributionsMultiNormalPrec, policiesMu) {
  boost::random::mt19937 rng;

  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();
  double ninf = -std::numeric_limits<double>::infinity();

  Eigen::VectorXd mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Eigen::MatrixXd Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;

  EXPECT_NO_THROW(stan::math::multi_normal_prec_rng(mu, Sigma, rng));

  mu(0) = inf;
  EXPECT_THROW(stan::math::multi_normal_prec_rng(mu, Sigma, rng),
               std::domain_error);

  mu(0) = ninf;
  EXPECT_THROW(stan::math::multi_normal_prec_rng(mu, Sigma, rng),
               std::domain_error);

  mu(0) = nan;
  EXPECT_THROW(stan::math::multi_normal_prec_rng(mu, Sigma, rng),
               std::domain_error);

  mu(0) = 1.0;
}

TEST(ProbDistributionsMultiNormalPrec, SizeMismatch) {
  boost::random::mt19937 rng;
  Eigen::VectorXd mu(2);
  mu << 1.0, -1.0;
  Eigen::MatrixXd Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  EXPECT_THROW(stan::math::multi_normal_prec_rng(mu, Sigma, rng),
               std::invalid_argument);
}

TEST(ProbDistributionsMultiNormalPrec, marginalOneChiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;

  Eigen::MatrixXd sigma(3, 3);
  sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 1.0, 0.0, 1.0, 3.0;
  Eigen::MatrixXd L = sigma.inverse();

  std::vector<Eigen::VectorXd> mu(3);
  mu[0].resize(3);
  mu[1].resize(3);
  mu[2].resize(3);
  mu[0] << 2.0, -2.0, 11.0;
  mu[1] << 7.0, -3.0, 5.0;
  mu[2] << 5.0, -6.0, 1.0;

  int N = 10000;
  int K = boost::math::round(2 * std::pow(N, 0.4));
  boost::math::normal_distribution<> dist(mu[0](0), sqrt(sigma(0, 0)));

  std::vector<double> quantiles;
  for (int i = 1; i < K; i++)
    quantiles.push_back(quantile(dist, i * std::pow(K, -1.0)));
  quantiles.push_back(std::numeric_limits<double>::max());

  Eigen::VectorXd a(mu[0].rows());
  std::vector<double> samples;
  for (int count = 0; count < N; ++count) {
    a = stan::math::multi_normal_prec_rng(mu, L, rng)[0];
    samples.push_back(a(0));
  }

  assert_matches_quantiles(samples, quantiles, 1e-6);
}

TEST(ProbDistributionsMultiNormalPrec, marginalTwoChiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;

  Eigen::MatrixXd sigma(3, 3);
  sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 1.0, 0.0, 1.0, 3.0;
  Eigen::MatrixXd L = sigma.inverse();

  std::vector<Eigen::RowVectorXd> mu(3);
  mu[0].resize(3);
  mu[1].resize(3);
  mu[2].resize(3);
  mu[0] << 1.0, 5.0, -1.0;
  mu[1] << 2.0, -2.0, 11.0;
  mu[2] << 7.0, 0.0, 3.0;

  int N = 10000;
  int K = boost::math::round(2 * std::pow(N, 0.4));
  boost::math::normal_distribution<> dist(mu[1](1), sqrt(sigma(1, 1)));

  std::vector<double> quantiles;
  for (int i = 1; i < K; i++)
    quantiles.push_back(quantile(dist, i * std::pow(K, -1.0)));
  quantiles.push_back(std::numeric_limits<double>::max());

  Eigen::VectorXd a(mu[0].rows());
  std::vector<double> samples;
  for (int count = 0; count < N; ++count) {
    a = stan::math::multi_normal_prec_rng(mu, L, rng)[1];
    samples.push_back(a(1));
  }

  assert_matches_quantiles(samples, quantiles, 1e-6);
}

TEST(ProbDistributionsMultiNormalPrec, marginalThreeChiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;

  Eigen::MatrixXd sigma(3, 3);
  sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 1.0, 0.0, 1.0, 16.0;
  Eigen::MatrixXd L = sigma.inverse();

  Eigen::VectorXd mu(3);
  mu << 2.0, -2.0, 11.0;

  int N = 10000;
  int K = boost::math::round(2 * std::pow(N, 0.4));
  boost::math::normal_distribution<> dist(mu(2), sqrt(sigma(2, 2)));

  std::vector<double> quantiles;
  for (int i = 1; i < K; i++)
    quantiles.push_back(quantile(dist, i * std::pow(K, -1.0)));
  quantiles.push_back(std::numeric_limits<double>::max());

  Eigen::VectorXd a(mu.rows());
  std::vector<double> samples;
  for (int count = 0; count < N; ++count) {
    a = stan::math::multi_normal_prec_rng(mu, L, rng);
    samples.push_back(a(2));
  }

  assert_matches_quantiles(samples, quantiles, 1e-6);
}
