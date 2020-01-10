#include <stan/math/prim.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/prob/VectorIntRNGTestRig.hpp>
#include <test/unit/math/prim/prob/util.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

class BernoulliTestRig : public VectorIntRNGTestRig {
 public:
  BernoulliTestRig()
      : VectorIntRNGTestRig(10000, 10, {0, 1}, {0.0, 0.1, 0.2, 0.7, 1.0},
                            {0, 1}, {-2.0, -0.5, 1.1, 2.0}, {-2, -1, 2}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& theta, const T2&, const T3&,
                        T_rng& rng) const {
    return stan::math::bernoulli_rng(theta, rng);
  }

  template <typename T1>
  double pmf(int y, T1 theta, double, double) const {
    return std::exp(stan::math::bernoulli_lpmf(y, theta));
  }
};

TEST(ProbDistributionsBernoulli, errorCheck) {
  check_dist_throws_all_types(BernoulliTestRig());
}

TEST(ProbDistributionsBernoulli, distributionCheck) {
  check_counts_real(BernoulliTestRig());
}

TEST(ProbDistributionsBernoulli, error_check) {
  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::bernoulli_rng(0.6, rng));

  EXPECT_THROW(stan::math::bernoulli_rng(1.6, rng), std::domain_error);
  EXPECT_THROW(stan::math::bernoulli_rng(-0.6, rng), std::domain_error);
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
  for (int i = 0; i < N; ++i) {
    ++counts[stan::math::bernoulli_rng(0.4, rng)];
  }

  assert_chi_squared(counts, expected, 1e-6);
}
