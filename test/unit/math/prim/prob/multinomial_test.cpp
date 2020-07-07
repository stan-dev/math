#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <limits>
#include <vector>

TEST(ProbDistributionsMultinomial, RNGSize) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> theta(5);
  // error in 2.1.0 due to overflow in binomial call due to division
  theta << 0.3, 0.1, 0.2, 0.2, 0.2;
  std::vector<int> sample = stan::math::multinomial_rng(theta, 10, rng);
  // bug in 2.1.0 returned 10 rather than 5 for returned size
  EXPECT_EQ(5U, sample.size());
}

TEST(ProbDistributionsMultinomial, Multinomial) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta << 0.2, 0.3, 0.5;
  EXPECT_FLOAT_EQ(-2.002481, stan::math::multinomial_log(ns, theta));
}
TEST(ProbDistributionsMultinomial, Propto) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta << 0.2, 0.3, 0.5;
  EXPECT_FLOAT_EQ(0.0, stan::math::multinomial_log<true>(ns, theta));
}

using stan::math::multinomial_log;

TEST(ProbDistributionsMultinomial, error) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();

  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta << 0.2, 0.3, 0.5;

  EXPECT_NO_THROW(multinomial_log(ns, theta));

  ns[1] = 0;
  EXPECT_NO_THROW(multinomial_log(ns, theta));
  ns[1] = -1;
  EXPECT_THROW(multinomial_log(ns, theta), std::domain_error);
  ns[1] = 1;

  theta(0) = 0.0;
  EXPECT_THROW(multinomial_log(ns, theta), std::domain_error);
  theta(0) = nan;
  EXPECT_THROW(multinomial_log(ns, theta), std::domain_error);
  theta(0) = inf;
  EXPECT_THROW(multinomial_log(ns, theta), std::domain_error);
  theta(0) = -inf;
  EXPECT_THROW(multinomial_log(ns, theta), std::domain_error);
  theta(0) = -1;
  theta(1) = 1.5;
  theta(2) = 0.5;
  EXPECT_THROW(multinomial_log(ns, theta), std::domain_error);
  theta(0) = 0.2;
  theta(1) = 0.3;
  theta(2) = 0.5;

  ns.resize(2);
  EXPECT_THROW(multinomial_log(ns, theta), std::invalid_argument);
}

TEST(ProbDistributionsMultinomial, zeros) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  double result;
  std::vector<int> ns;
  ns.push_back(0);
  ns.push_back(1);
  ns.push_back(2);
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta << 0.2, 0.3, 0.5;

  result = multinomial_log(ns, theta);
  EXPECT_FALSE(std::isnan(result));

  std::vector<int> ns2;
  ns2.push_back(0);
  ns2.push_back(0);
  ns2.push_back(0);

  double result2 = multinomial_log(ns2, theta);
  EXPECT_FLOAT_EQ(0.0, result2);
}

TEST(ProbDistributionsMultinomial, error_check) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  boost::random::mt19937 rng;

  Matrix<double, Dynamic, 1> theta(3);
  theta << 0.15, 0.45, 0.40;

  EXPECT_THROW(stan::math::multinomial_rng(theta, -3, rng), std::domain_error);

  theta << 0.15, 0.45, 0.50;
  EXPECT_THROW(stan::math::multinomial_rng(theta, 3, rng), std::domain_error);
}

TEST(ProbDistributionsMultinomial, chiSquareGoodnessFitTest) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  boost::random::mt19937 rng;
  int M = 10;
  int trials = 1000;
  int N = M * trials;

  int K = 3;
  Matrix<double, Dynamic, 1> theta(K);
  theta << 0.2, 0.35, 0.45;
  boost::math::chi_squared mydist(K - 1);

  double expect[K];
  for (int i = 0; i < K; ++i)
    expect[i] = N * theta(i);

  int bin[K];
  for (int i = 0; i < K; ++i)
    bin[i] = 0;

  for (int count = 0; count < M; ++count) {
    std::vector<int> a = stan::math::multinomial_rng(theta, trials, rng);
    for (int i = 0; i < K; ++i)
      bin[i] += a[i];
  }

  double chi = 0;
  for (int j = 0; j < K; j++)
    chi += ((bin[j] - expect[j]) * (bin[j] - expect[j])) / expect[j];

  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}
