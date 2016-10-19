#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <test/unit/math/prim/scal/prob/util.hpp>

TEST(ProbDistributionsBernoulli, error_check) {
  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::bernoulli_rng(0.6, rng));

  EXPECT_THROW(stan::math::bernoulli_rng(1.6, rng),std::domain_error);
  EXPECT_THROW(stan::math::bernoulli_rng(-0.6, rng),std::domain_error);
  EXPECT_THROW(stan::math::bernoulli_rng(stan::math::positive_infinity(), rng),
               std::domain_error);
}

TEST(ProbDistributionsBernoulli, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000;

  std::vector<double> expected;
  expected.push_back(N * (1 - 0.4));
  expected.push_back(N * 0.4);

  std::vector<int> counts(2);
  for (int i=0; i<N; ++i) {
    ++counts[stan::math::bernoulli_rng(0.4, rng)];
  }

  assert_chi_squared(counts, expected, 1e-6);
}
