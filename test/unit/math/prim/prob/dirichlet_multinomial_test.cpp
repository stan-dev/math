#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <limits>
#include <vector>

using stan::math::dirichlet_multinomial_log;
using stan::math::dirichlet_multinomial_rng;

TEST(ProbDistributionsDirichletMultinomial, DirichletMultinomial) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta << 2.0, 3.0, 5.0;
  EXPECT_FLOAT_EQ(-2.477938, dirichlet_multinomial_log(ns, theta));
}

TEST(ProbDistributionsDirichletMultinomial, Propto) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta << 2.0, 3.0, 5.0;
  EXPECT_FLOAT_EQ(0.0, dirichlet_multinomial_log<true>(ns, theta));
}

TEST(ProbDistributionsDirichletMultinomial, error) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();

  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta << 2.0, 3.0, 5.0;

  EXPECT_NO_THROW(dirichlet_multinomial_log(ns, theta));

  ns[1] = 0;
  EXPECT_NO_THROW(dirichlet_multinomial_log(ns, theta));
  ns[1] = -1;
  EXPECT_THROW(dirichlet_multinomial_log(ns, theta), std::domain_error);
  ns[1] = 1;

  theta(0) = 0.0;
  EXPECT_THROW(dirichlet_multinomial_log(ns, theta), std::domain_error);
  theta(0) = -1.0;
  EXPECT_THROW(dirichlet_multinomial_log(ns, theta), std::domain_error);
  theta(0) = nan;
  EXPECT_THROW(dirichlet_multinomial_log(ns, theta), std::domain_error);
  theta(0) = inf;
  EXPECT_THROW(dirichlet_multinomial_log(ns, theta), std::domain_error);
  theta(0) = -inf;
  EXPECT_THROW(dirichlet_multinomial_log(ns, theta), std::domain_error);
  theta(0) = 2.0;

  ns.resize(2);
  EXPECT_THROW(dirichlet_multinomial_log(ns, theta), std::invalid_argument);
}

TEST(ProbDistributionsDirichletMultinomial, zeros) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  std::vector<int> ns;
  ns.push_back(0);
  ns.push_back(0);
  ns.push_back(0);
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta << 2.0, 3.0, 5.0;

  double result = dirichlet_multinomial_log(ns, theta);
  EXPECT_FLOAT_EQ(0.0, result);
}

TEST(ProbDistributionsDirichletMultinomial, rng) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  boost::random::mt19937 rng;

  Matrix<double, Dynamic, 1> theta(3);
  theta << 0.5, 5.5, 12.0;

  int n = 100;

  std::vector<int> a = dirichlet_multinomial_rng(theta, n, rng);

  EXPECT_EQ(std::accumulate(a.begin(), a.end(), 0), n);
  for (int k : a) EXPECT_GE(k, 0);

  std::vector<int> b = dirichlet_multinomial_rng(theta, 0, rng);
  for (int k : b) EXPECT_EQ(k, 0);
}

TEST(ProbDistributionsDirichletMultinomial, rngError) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  boost::random::mt19937 rng;

  Matrix<double, Dynamic, 1> theta(3);
  theta << 0.5, 1.5, 4.0;

  EXPECT_THROW(dirichlet_multinomial_rng(theta, -3, rng), std::domain_error);

  EXPECT_NO_THROW(dirichlet_multinomial_rng(theta, 0, rng));

  theta << 0.0, 2.5, 5.0;
  EXPECT_THROW(dirichlet_multinomial_rng(theta, 3, rng), std::domain_error);

  theta << 1.2, -2.8, 5.0;
  EXPECT_THROW(dirichlet_multinomial_rng(theta, 5, rng), std::domain_error);
}

TEST(ProbDistributionsDirichletMultinomial, chiSquareGoodnessFitTest) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  boost::random::mt19937 rng;
  int M = 100;
  int trials = 100;
  int N = M * trials;

  int K = 4;
  Matrix<double, Dynamic, 1> theta(K);
  theta << 2.0, 0.5, 4.5, 10.0;
  boost::math::chi_squared mydist(K - 1);

  double expect[K];
  for (int i = 0; i < K; ++i)
    expect[i] = N * theta(i) / theta.sum();

  int bin[K];
  for (int i = 0; i < K; ++i)
    bin[i] = 0;

  for (int count = 0; count < M; ++count) {
    std::vector<int> a = dirichlet_multinomial_rng(theta, trials, rng);
    for (int i = 0; i < K; ++i)
      bin[i] += a[i];
  }

  double chi = 0;
  for (int j = 0; j < K; j++)
    chi += ((bin[j] - expect[j]) * (bin[j] - expect[j])) / expect[j];

  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}
