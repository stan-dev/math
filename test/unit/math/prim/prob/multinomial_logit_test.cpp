#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <limits>
#include <vector>

using Eigen::Dynamic;
using Eigen::Matrix;

TEST(ProbDistributionsMultinomialLogit, RNGSize) {
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> beta(5);
  beta << log(0.3), log(0.1), log(0.2), log(0.2), log(0.2);
  std::vector<int> sample = stan::math::multinomial_logit_rng(beta, 10, rng);
  EXPECT_EQ(5U, sample.size());
}

TEST(ProbDistributionsMultinomialLogit, MultinomialLogit) {
  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Matrix<double, Dynamic, 1> beta(3, 1);
  beta << log(0.2), log(0.3), log(0.5);
  EXPECT_FLOAT_EQ(-2.002481, stan::math::multinomial_logit_log(ns, beta));
}
TEST(ProbDistributionsMultinomialLogit, Propto) {
  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Matrix<double, Dynamic, 1> beta(3, 1);
  beta << log(0.2), log(0.3), log(0.5);
  EXPECT_FLOAT_EQ(0.0, stan::math::multinomial_logit_log<true>(ns, beta));
}

using stan::math::multinomial_logit_log;

TEST(ProbDistributionsMultinomialLogit, error) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();

  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Matrix<double, Dynamic, 1> beta(3, 1);
  beta << log(0.2), log(0.3), log(0.5);

  EXPECT_NO_THROW(multinomial_logit_log(ns, beta));

  ns[1] = 0;
  EXPECT_NO_THROW(multinomial_logit_log(ns, beta));
  ns[1] = -1;
  EXPECT_THROW(multinomial_logit_log(ns, beta), std::domain_error);
  ns[1] = 1;

  beta(0) = nan;
  EXPECT_THROW(multinomial_logit_log(ns, beta), std::domain_error);
  beta(0) = inf;
  EXPECT_THROW(multinomial_logit_log(ns, beta), std::domain_error);
  beta(0) = -inf;
  EXPECT_THROW(multinomial_logit_log(ns, beta), std::domain_error);

  beta(0) = 0.2;
  beta(1) = 0.3;
  beta(2) = 0.5;

  ns.resize(2);
  EXPECT_THROW(multinomial_logit_log(ns, beta), std::invalid_argument);
}

TEST(ProbDistributionsMultinomialLogit, zeros) {
  double result;
  std::vector<int> ns;
  ns.push_back(0);
  ns.push_back(1);
  ns.push_back(2);
  Matrix<double, Dynamic, 1> beta(3, 1);
  beta << log(0.2), log(0.3), log(0.5);

  result = multinomial_logit_log(ns, beta);
  EXPECT_FALSE(std::isnan(result));

  std::vector<int> ns2;
  ns2.push_back(0);
  ns2.push_back(0);
  ns2.push_back(0);

  double result2 = multinomial_logit_log(ns2, beta);
  EXPECT_FLOAT_EQ(0.0, result2);
}

TEST(ProbDistributionsMultinomialLogit, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int M = 10;
  int trials = 1000;
  int N = M * trials;

  int K = 3;
  Matrix<double, Dynamic, 1> beta(K);
  beta << -1, 1, -10;
  Eigen::VectorXd theta = stan::math::softmax(beta);
  boost::math::chi_squared mydist(K - 1);

  double expect[K];
  for (int i = 0; i < K; ++i)
    expect[i] = N * theta(i);

  int bin[K];
  for (int i = 0; i < K; ++i)
    bin[i] = 0;

  for (int count = 0; count < M; ++count) {
    std::vector<int> a = stan::math::multinomial_logit_rng(beta, trials, rng);
    for (int i = 0; i < K; ++i)
      bin[i] += a[i];
  }

  double chi = 0;
  for (int j = 0; j < K; j++)
    chi += ((bin[j] - expect[j]) * (bin[j] - expect[j])) / expect[j];

  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}
