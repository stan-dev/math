#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <test/unit/math/prim/scal/prob/util.hpp>

TEST(ProbDistributionsBernoulliLogit, error_check) {
  boost::random::mt19937 rng;

  EXPECT_NO_THROW(stan::math::bernoulli_logit_rng(-3.5, rng));
  EXPECT_THROW(stan::math::bernoulli_logit_rng(stan::math::positive_infinity(), rng),
               std::domain_error);
}

TEST(ProbDistributionsBernoulliLogit, logitChiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000; // number of samples
 
  double parameter = -0.5; // logit-transformed probability
  double prob = stan::math::inv_logit(-0.5); // actual probability

  std::vector<double> expected;
  expected.push_back(N * (1 - prob));
  expected.push_back(N * prob);

  std::vector<int> counts(2);
  for (int i=0; i<N; ++i) {
    ++counts[stan::math::bernoulli_logit_rng(parameter, rng)];
  }

  assert_chi_squared(counts, expected, 1e-6);
}
