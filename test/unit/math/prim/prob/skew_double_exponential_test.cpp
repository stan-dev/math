#include <stan/math/prim.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/prob/util.hpp>
#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <limits>
#include <vector>

class SkewDoubleExponentialTestRig : public VectorRNGTestRig {
 public:
  SkewDoubleExponentialTestRig()
      : VectorRNGTestRig(10000, 10, {-2.5, -1.7, -0.1, 0.1, 2.0},
                         {-3, -2, -1, 0, 2, 6}, {}, {}, {0.1, 1.0, 2.5, 4.0},
                         {1, 2, 3, 4}, {-2.7, -1.5, -0.5, 0.0}, {-3, -2, -1, 0},
                         {0.5}, {0}, {-0.1, 1.1}, {}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& mu, const T2& sigma, const T3& tau,
                        T_rng& rng) const {
    return stan::math::skew_double_exponential_rng(mu, sigma, tau, rng);
  }
};

double icdf(double z, double mu, double sigma, double tau) {
  if (z < tau) {
    return log(z / tau) * sigma / (2.0 * (1.0 - tau)) + mu;
  } else {
    return log((1.0 - z) / (1.0 - tau)) * (-sigma) / (2.0 * tau) + mu;
  }
}

TEST(ProbDistributionsSkewedDoubleExponential, errorCheck) {
  check_dist_throws_all_types(SkewDoubleExponentialTestRig());
}

TEST(ProbDistributionsSkewedDoubleExponential, error_check) {
  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::skew_double_exponential_rng(10.0, 2.0, .1, rng));

  EXPECT_THROW(stan::math::skew_double_exponential_rng(10.0, 2.0, -1.0, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::skew_double_exponential_rng(
                   10, 2, stan::math::positive_infinity(), rng),
               std::domain_error);
}

TEST(ProbDistributionsSkewedDoubleExponential, test_sampling_icdf) {
  for (double p : {0.0, 0.1, 0.2, 0.5, 0.7, 0.9, 0.99}) {
    for (double mu : {-1.11, 0.13, 1.2, 4.67}) {
      for (double sigma : {0.11, 1.33}) {
        for (double tau : {0.1, 0.4, 0.77, 0.89}) {
          double x = icdf(p, mu, sigma, tau);
          EXPECT_FLOAT_EQ(
              stan::math::skew_double_exponential_cdf(x, mu, sigma, tau), p);
        }
      }
    }
  }
}

TEST(ProbDistributionsSkewedDoubleExponential, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000;
  int K = stan::math::round(2 * std::pow(N, 0.4));

  std::vector<double> samples;
  for (int i = 0; i < N; ++i) {
    samples.push_back(
        stan::math::skew_double_exponential_rng(2.0, 1.0, 0.25, rng));
  }
  std::vector<double> quantiles;
  for (int i = 1; i < K; ++i) {
    double frac = static_cast<double>(i) / K;
    quantiles.push_back(icdf(frac, 2.0, 1.0, 0.25));
  }
  quantiles.push_back(std::numeric_limits<double>::max());

  // Assert that they match
  assert_matches_quantiles(samples, quantiles, 1e-6);
}
