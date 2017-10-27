#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim/scal.hpp>

TEST(ProbDistributionBetaBinomial, error_check) {
  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::beta_binomial_rng(4, 0.6, 2.0, rng));

  EXPECT_THROW(stan::math::beta_binomial_rng(-4, 0.6, 2, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::beta_binomial_rng(4, -0.6, 2, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::beta_binomial_rng(4, 0.6, -2, rng),
               std::domain_error);
  EXPECT_THROW(
      stan::math::beta_binomial_rng(4, stan::math::positive_infinity(), 2, rng),
      std::domain_error);
  EXPECT_THROW(stan::math::beta_binomial_rng(
                   4, 0.6, stan::math::positive_infinity(), rng),
               std::domain_error);
}
