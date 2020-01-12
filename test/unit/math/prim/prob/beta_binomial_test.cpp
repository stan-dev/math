#include <stan/math/prim.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/prob/VectorIntRNGTestRig.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

class BetaBinomialTestRig : public VectorIntRNGTestRig {
 public:
  BetaBinomialTestRig()
      : VectorIntRNGTestRig(10000, 10, {0, 1, 2, 3, 4, 5, 6}, {}, {0, 1, 3, 8},
                            {}, {-1, -5, -7}, {0.1, 1.7, 3.99}, {1, 2, 3},
                            {-2.1, -0.5, 0.0}, {-3, -1, 0}, {0.1, 1.1, 4.99},
                            {1, 2, 3}, {-3.0, -2.0, 0.0}, {-3, -1, 0}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& N, const T2& alpha, const T3& beta,
                        T_rng& rng) const {
    return stan::math::beta_binomial_rng(N, alpha, beta, rng);
  }

  template <typename T1>
  double pmf(int y, T1 N, double alpha, double beta) const {
    if (y <= N) {
      return std::exp(stan::math::beta_binomial_lpmf(y, N, alpha, beta));
    } else {
      return 0.0;
    }
  }
};

TEST(ProbDistributionsBetaBinomial, errorCheck) {
  check_dist_throws_int_first_argument(BetaBinomialTestRig());
}

TEST(ProbDistributionsBetaBinomial, distributionCheck) {
  check_counts_int_real_real(BetaBinomialTestRig());
}

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
