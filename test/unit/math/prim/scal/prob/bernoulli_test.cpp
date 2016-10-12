#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>

TEST(ProbDistributionsBernoulli, error_check) {
  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::bernoulli_rng(0.6,rng));
  EXPECT_NO_THROW(stan::math::bernoulli_logit_rng(-3.5,rng));

  EXPECT_THROW(stan::math::bernoulli_rng(1.6,rng),std::domain_error);
  EXPECT_THROW(stan::math::bernoulli_rng(-0.6,rng),std::domain_error);
  EXPECT_THROW(stan::math::bernoulli_rng(stan::math::positive_infinity(),rng),
               std::domain_error);
  EXPECT_THROW(stan::math::bernoulli_logit_rng(stan::math::positive_infinity(),rng),
               std::domain_error);
}

/*
 * Uses a chi-squared test to check that counts are reasonably consistent
 * with draws from a Bernoulli distribution..
 */
void assert_bernoulli(const std::vector<int>& counts, const std::vector<double>& expected_counts, double tolerance) {
  boost::math::chi_squared dist(1);
  
  EXPECT_EQ(counts.size(), 2);
  EXPECT_EQ(expected_counts.size(), 2);
  
  double chi = 0;
  for (int i=0; i<2; ++i) {
    const double discrepancy = expected_counts[i] - counts[i];
    chi += discrepancy * discrepancy / expected_counts[i];
  }
  
  const double chi_threshold = quantile(complement(dist, tolerance));
  EXPECT_TRUE(chi < chi_threshold);
}

TEST(ProbDistributionsBernoulli, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000;

  const int tmp1[2] = {0, 0};
  std::vector<int> counts(tmp1, tmp1 + 2);

  const double tmp2[2] = {N * (1 - 0.4), N * 0.4};
  const std::vector<double> expected_counts(tmp2, tmp2 + 2);

  for (int i=0; i<N; ++i) {
    ++counts[stan::math::bernoulli_rng(0.4, rng)];
  }
  assert_bernoulli(counts, expected_counts, 1e-6);
}


TEST(ProbDistributionsBernoulli, logitChiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000; // number of samples
 
  double parameter = -0.5; // logit-transformed probability
  double prob = stan::math::inv_logit(-0.5); // actual probability

  const int tmp1[2] = {0, 0};
  std::vector<int> counts(tmp1, tmp1 + 2);

  const double tmp2[2] = {N * (1 - prob), N * prob};
  const std::vector<double> expected_counts(tmp2, tmp2 + 2);

  for (int i=0; i<N; ++i) {
    ++counts[stan::math::bernoulli_logit_rng(parameter, rng)];
  }
  assert_bernoulli(counts, expected_counts, 1e-6);
}
