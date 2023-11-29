#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <limits>
#include <vector>

using stan::math::dirichlet_multinomial_lpmf;
using stan::math::dirichlet_multinomial_rng;

TEST(ProbDistributionsDirichletMultinomial, DirichletMultinomial) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  std::vector<int> ns = {1, 2, 3};
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta << 2.0, 3.0, 5.0;
  // test against some pre-computed log-prob values
  EXPECT_FLOAT_EQ(-2.477938, dirichlet_multinomial_lpmf(ns, theta));

  ns = {2, 12, 31};
  theta << 0.2, 6.1, 3.0;
  EXPECT_FLOAT_EQ(-8.1557148, dirichlet_multinomial_lpmf(ns, theta));

  ns = {106, 1203, 31020};
  theta << 50.0, 40.0, 20.0;
  EXPECT_FLOAT_EQ(-304.5664433, dirichlet_multinomial_lpmf(ns, theta));
}

TEST(ProbDistributionsDirichletMultinomial, Propto) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  std::vector<int> ns = {1, 2, 3};
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta << 2.0, 3.0, 5.0;
  // if propto is true, then result should be 0.0
  EXPECT_FLOAT_EQ(0.0, dirichlet_multinomial_lpmf<true>(ns, theta));
}

TEST(ProbDistributionsDirichletMultinomial, error) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();

  std::vector<int> ns = {1, 2, 3};
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta << 2.0, 3.0, 5.0;

  EXPECT_NO_THROW(dirichlet_multinomial_lpmf(ns, theta));

  // zero ns are fine
  ns[1] = 0;
  EXPECT_NO_THROW(dirichlet_multinomial_lpmf(ns, theta));

  // negative values for ns are not allowed
  ns[1] = -1;
  EXPECT_THROW(dirichlet_multinomial_lpmf(ns, theta), std::domain_error);
  ns[1] = 1;

  // prior sizes should be strictly positive, finite, and not NaN
  theta(0) = 0.0;
  EXPECT_THROW(dirichlet_multinomial_lpmf(ns, theta), std::domain_error);
  theta(0) = -1.0;
  EXPECT_THROW(dirichlet_multinomial_lpmf(ns, theta), std::domain_error);
  theta(0) = nan;
  EXPECT_THROW(dirichlet_multinomial_lpmf(ns, theta), std::domain_error);
  theta(0) = inf;
  EXPECT_THROW(dirichlet_multinomial_lpmf(ns, theta), std::domain_error);
  theta(0) = -inf;
  EXPECT_THROW(dirichlet_multinomial_lpmf(ns, theta), std::domain_error);
  theta(0) = 2.0;

  // the size of ns and theta should be identical
  ns.resize(2);
  EXPECT_THROW(dirichlet_multinomial_lpmf(ns, theta), std::invalid_argument);
}

TEST(ProbDistributionsDirichletMultinomial, zeros) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  std::vector<int> ns(3, 0);
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta << 2.0, 3.0, 5.0;

  // if all ns are zero, then the log prob is zero
  double result = dirichlet_multinomial_lpmf(ns, theta);
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

  // PRNG samples should sum to n, and each element is non-negative
  EXPECT_EQ(std::accumulate(a.begin(), a.end(), 0), n);
  for (int k : a)
    EXPECT_GE(k, 0);

  // if n is zero, then each element of the sample is zero
  std::vector<int> b = dirichlet_multinomial_rng(theta, 0, rng);
  for (int k : b)
    EXPECT_EQ(k, 0);
}

TEST(ProbDistributionsDirichletMultinomial, rngError) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  boost::random::mt19937 rng;

  Matrix<double, Dynamic, 1> theta(3);
  theta << 0.5, 1.5, 4.0;

  // number of trials parameter should be non-negative
  EXPECT_THROW(dirichlet_multinomial_rng(theta, -3, rng), std::domain_error);

  // but zero trials is a valid argument
  EXPECT_NO_THROW(dirichlet_multinomial_rng(theta, 0, rng));

  // prior sizes should be strictly positive
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

TEST(ProbDistributionsDirichletMultinomial, equivBetaBinomial) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::beta_binomial_lpmf;

  boost::random::mt19937 rng;
  std::vector<int> ns;
  Matrix<double, Dynamic, 1> theta(2, 1);
  double alpha = 4.0;
  double beta = 5.0;
  theta << alpha, beta;

  // number of random samples from DirMult distribution
  int draws = 100;
  // population size parameter
  int N = 25;

  for (int i = 0; i < draws; ++i) {
    ns = dirichlet_multinomial_rng(theta, N, rng);
    // For K = 2, BetaBinom and DirMult should be identical
    EXPECT_FLOAT_EQ(beta_binomial_lpmf(ns[0], N, alpha, beta),
                    dirichlet_multinomial_lpmf(ns, theta));
  }

  // test again with small prior sizes
  alpha = 0.5;
  beta = 0.1;
  theta << alpha, beta;

  for (int i = 0; i < draws; ++i) {
    ns = dirichlet_multinomial_rng(theta, N, rng);
    // For K = 2, BetaBinom and DirMult should be identical
    EXPECT_FLOAT_EQ(beta_binomial_lpmf(ns[0], N, alpha, beta),
                    dirichlet_multinomial_lpmf(ns, theta));
  }
}
