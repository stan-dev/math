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

TEST(ProbDistributionsBernoulli, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000;
  boost::math::bernoulli_distribution<>dist (0.4);
  boost::math::chi_squared mydist(1);
 
  int bin[2] = {0, 0};
  double expect [2] = {N * (1 - 0.4), N * (0.4)};

  int count = 0;

  while (count < N) {
    int a = stan::math::bernoulli_rng(0.4,rng);
    if(a == 1)
      ++bin[1];
    else
      ++bin[0];
    count++;
   }

  double chi = 0;

  for(int j = 0; j < 2; j++)
    chi += ((bin[j] - expect[j]) * (bin[j] - expect[j]) / expect[j]);

  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}


TEST(ProbDistributionsBernoulli, logitChiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000;
  boost::math::bernoulli_distribution<>dist (0.4);
  boost::math::chi_squared mydist(1);
 
  double parameter = -0.5;
  double prob = stan::math::inv_logit(-0.5);
  int bin[2] = {0, 0};
  double expect [2] = {N * (1 - prob), N * prob};

  int count = 0;

  while (count < N) {
    int a = stan::math::bernoulli_logit_rng(parameter,rng);
    if(a == 1)
      ++bin[1];
    else
      ++bin[0];
    count++;
   }

  double chi = 0;

  for(int j = 0; j < 2; j++)
    chi += ((bin[j] - expect[j]) * (bin[j] - expect[j]) / expect[j]);

  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}
