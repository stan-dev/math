#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <limits>
#include <vector>

using stan::math::poisson_binomial_lccdf;
using stan::math::poisson_binomial_lcdf;
using stan::math::poisson_binomial_lpmf;
using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

static double inff = std::numeric_limits<double>::infinity();

TEST(ProbDistributionsPoissonBinomial,
     poisson_binomial_check_error_scalar_y_oob) {
  vec theta(3);
  theta << 0.5, 0.2, 0.7;

  EXPECT_THROW(poisson_binomial_lpmf(-1, theta), std::domain_error);
  EXPECT_THROW(poisson_binomial_lpmf(4, theta), std::domain_error);
  EXPECT_THROW(poisson_binomial_lcdf(-1, theta), std::domain_error);
  EXPECT_THROW(poisson_binomial_lcdf(4, theta), std::domain_error);
  EXPECT_THROW(poisson_binomial_lccdf(-1, theta), std::domain_error);
  EXPECT_THROW(poisson_binomial_lccdf(4, theta), std::domain_error);
}

TEST(ProbDistributionsPoissonBinomial, check_error_theta_is_not_prob) {
  vec theta(3);
  theta << 0.5, 0.2, 1.1;
  vec theta2(3);
  theta2 << 0.5, 0.2, inff;

  EXPECT_THROW(poisson_binomial_lpmf(1, theta), std::domain_error);
  EXPECT_THROW(poisson_binomial_lpmf(1, theta2), std::domain_error);
  EXPECT_THROW(poisson_binomial_lcdf(1, theta), std::domain_error);
  EXPECT_THROW(poisson_binomial_lcdf(1, theta2), std::domain_error);
  EXPECT_THROW(poisson_binomial_lccdf(1, theta), std::domain_error);
  EXPECT_THROW(poisson_binomial_lccdf(1, theta2), std::domain_error);
}

TEST(ProbDistributionsPoissonBinomial, check_error_vectorial_y_oob) {
  vec theta(3);
  theta << 0.5, 0.2, 0.1;
  std::vector<int> ys1{-1, 2};
  std::vector<int> ys2{4, 3};

  EXPECT_THROW(poisson_binomial_lpmf(ys1, theta), std::domain_error);
  EXPECT_THROW(poisson_binomial_lpmf(ys2, theta), std::domain_error);
  EXPECT_THROW(poisson_binomial_lcdf(ys1, theta), std::domain_error);
  EXPECT_THROW(poisson_binomial_lcdf(ys2, theta), std::domain_error);
  EXPECT_THROW(poisson_binomial_lccdf(ys1, theta), std::domain_error);
  EXPECT_THROW(poisson_binomial_lccdf(ys2, theta), std::domain_error);
}

TEST(ProbDistributionsPoissonBinomial,
     check_error_vectorial_y_theta_is_not_prob) {
  vec theta(3);
  theta << inff, 0.2, 0.1;
  std::vector<int> y{0, 2};

  EXPECT_THROW(poisson_binomial_lpmf(y, theta), std::domain_error);
  EXPECT_THROW(poisson_binomial_lcdf(y, theta), std::domain_error);
  EXPECT_THROW(poisson_binomial_lccdf(y, theta), std::domain_error);
}

TEST(ProbDistributionsPoissonBinomial,
     check_error_vectorial_theta_is_not_prob) {
  vec theta1(3);
  theta1 << -0.1, 0.2, 0.1;
  vec theta2(3);
  theta2 << 0.5, 0.2, 1.1;
  vec theta3(3);
  theta3 << 0.5, 0.2, inff;

  std::vector<int> ys{0, 2, 3};
  std::vector<vec> thetas{theta1, theta2, theta3};

  EXPECT_THROW(poisson_binomial_lpmf(ys, thetas), std::domain_error);
  EXPECT_THROW(poisson_binomial_lcdf(ys, thetas), std::domain_error);
  EXPECT_THROW(poisson_binomial_lccdf(ys, thetas), std::domain_error);
}

TEST(ProbDistributionsPoissonBinomial,
     check_error_vectorial_y_oob_with_vectorial_theta) {
  vec theta1(3);
  theta1 << 0.5, 0.2, 0.1;
  vec theta2(3);
  theta2 << 0.5, 0.2, 0.1;

  std::vector<int> ys{-1, 4};
  std::vector<vec> thetas{theta1, theta2};

  EXPECT_THROW(poisson_binomial_lpmf(ys, thetas), std::domain_error);
  EXPECT_THROW(poisson_binomial_lcdf(ys, thetas), std::domain_error);
  EXPECT_THROW(poisson_binomial_lccdf(ys, thetas), std::domain_error);
}

TEST(ProbDistributionsPoissonBinomial,
     check_error_vectorial_y_and_theta_inconsistent_sizes) {
  vec theta1(3);
  theta1 << 0.5, 0.2, 0.1;
  vec theta2(3);
  theta2 << 0.5, 0.2, 0.1;

  std::vector<int> ys{0};
  std::vector<vec> thetas{theta1, theta2};

  EXPECT_THROW(poisson_binomial_lpmf(ys, thetas), std::invalid_argument);
  EXPECT_THROW(poisson_binomial_lcdf(ys, thetas), std::invalid_argument);
  EXPECT_THROW(poisson_binomial_lccdf(ys, thetas), std::invalid_argument);
}

/*
 * Since there is no easy way to compute Poisson binomial quantiles with Boost,
 * test the distribution by emulating a binomial, i.e., by setting all
 * success probabilities to the same value
 */
TEST(ProbDistributionsPoissonBinomial, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000;
  int K = stan::math::round(2 * std::pow(N, 0.4));
  boost::math::binomial_distribution<> dist(100, 0.6);
  boost::math::chi_squared mydist(K - 1);

  int loc[K - 1];
  for (int i = 1; i < K; i++) {
    loc[i - 1] = i - 1;
  }

  int count = 0;
  int bin[K];
  double expect[K];
  for (int i = 0; i < K; i++) {
    bin[i] = 0;
    expect[i] = N * pdf(dist, i);
  }
  expect[K - 1] = N * (1 - cdf(dist, K - 2));

  vec probs(100);
  probs.fill(0.6);
  while (count < N) {
    int a = stan::math::poisson_binomial_rng(probs, rng);
    int i = 0;
    while (i < K - 1 && a > loc[i]) {
      ++i;
    }
    ++bin[i];
    count++;
  }

  double chi = 0;
  for (int j = 0; j < K; j++) {
    chi += ((bin[j] - expect[j]) * (bin[j] - expect[j]) / expect[j]);
  }

  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}
