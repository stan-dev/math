#include <stan/math/prim.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/prob/util.hpp>
#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <limits>
#include <vector>

class LoglogisticTestRig : public VectorRNGTestRig {
 public:
  LoglogisticTestRig()
      : VectorRNGTestRig(10000, 10, {2.5, 1.7, 0.2, 0.1, 2.0},
                         {3, 2, 1, 5, 10, 6}, {-2.5, -1.7, -0.2, -0.1, 0.0},
                         {-3, -2, -1, -4, -10, 0}, {0.1, 1.0, 2.5, 4.0},
                         {1, 2, 3, 4}, {-2.7, -1.5, -0.5, 0.0},
                         {-3, -2, -1, 0}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& alpha, const T2& beta, const T3& unused,
                        T_rng& rng) const {
    return stan::math::loglogistic_rng(alpha, beta, rng);
  }
};

double icdf(double x, double alpha, double beta) {
  return alpha * pow(x / (1 - x), 1 / beta);
}

TEST(ProbDistributionsLoglogistic, errorCheck) {
  check_dist_throws_all_types(LoglogisticTestRig());
}

TEST(ProbDistributionsLoglogistic, error_check) {
  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::loglogistic_rng(10.0, 2.0, rng));

  EXPECT_THROW(stan::math::loglogistic_rng(2.0, -1.0, rng), std::domain_error);
  EXPECT_THROW(stan::math::loglogistic_rng(-2.0, 1.0, rng), std::domain_error);
  EXPECT_THROW(
      stan::math::loglogistic_rng(10, stan::math::positive_infinity(), rng),
      std::domain_error);
  EXPECT_THROW(
      stan::math::loglogistic_rng(stan::math::positive_infinity(), 2, rng),
      std::domain_error);
}

TEST(ProbDistributionsLoglogistic, test_sampling_icdf) {
  for (double p : {0.0, 0.1, 0.2, 0.5, 0.7, 0.9, 0.99}) {
    for (double alpha : {1.11, 0.13, 1.2, 4.67}) {
      for (double beta : {0.11, 1.33, 2.0, 3.2}) {
        double x = icdf(p, alpha, beta);
        EXPECT_FLOAT_EQ(stan::math::loglogistic_cdf(x, alpha, beta), p);
      }
    }
  }
}

TEST(ProbDistributionsLoglogistic, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000;
  int K = stan::math::round(2 * std::pow(N, 0.4));

  std::vector<double> samples;
  for (int i = 0; i < N; ++i) {
    samples.push_back(stan::math::loglogistic_rng(1.2, 2.0, rng));
  }
  std::vector<double> quantiles;
  for (int i = 1; i < K; ++i) {
    double frac = static_cast<double>(i) / K;
    quantiles.push_back(icdf(frac, 1.2, 2.0));
  }
  quantiles.push_back(std::numeric_limits<double>::max());

  // Assert that they match
  assert_matches_quantiles(samples, quantiles, 1e-6);
}
