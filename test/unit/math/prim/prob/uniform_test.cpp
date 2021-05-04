#include <stan/math/prim.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/prob/util.hpp>
#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <limits>
#include <vector>

class UniformTestRig : public VectorRealRNGTestRig {
 public:
  UniformTestRig(std::vector<double> good_p1, std::vector<int> good_p1_int,
                 std::vector<double> bad_p1, std::vector<int> bad_p1_int,
                 std::vector<double> good_p2, std::vector<int> good_p2_int,
                 std::vector<double> bad_p2, std::vector<int> bad_p2_int)
      : VectorRealRNGTestRig(10000, 10, good_p1, good_p1_int, bad_p1,
                             bad_p1_int, good_p2, good_p2_int, bad_p2,
                             bad_p2_int) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& alpha, const T2& beta, const T3& unused,
                        T_rng& rng) const {
    return stan::math::uniform_rng(alpha, beta, rng);
  }

  std::vector<double> generate_quantiles(double alpha, double beta,
                                         double unused) const {
    std::vector<double> quantiles;
    double K = stan::math::round(2 * std::pow(N_, 0.4));
    boost::math::uniform_distribution<> dist(alpha, beta);

    for (int i = 1; i < K; ++i) {
      double frac = i / K;
      quantiles.push_back(quantile(dist, frac));
    }
    quantiles.push_back(std::numeric_limits<double>::max());

    return quantiles;
  }
};

TEST(ProbDistributionsUniform, errorCheck) {
  // beta must be greater than alpha, so we have to be careful about how the
  // tests are initialized here (the VectorRNGTestRig class assumes no order)

  check_dist_throws_all_types(
      UniformTestRig({-2.5, -1.7, -0.1, 0.0, 0.5}, {-3, -2, -1, 0, 1}, {}, {},
                     {1.1, 2.2, 3.8}, {2, 3, 6}, {}, {}));

  check_dist_throws_all_types(
      UniformTestRig({-7.5, -6.7, -5.1, -4.0}, {-7, -6, -5, -4}, {}, {},
                     {-3.0, -2.2, 0.0, 1.0}, {-3, -2, 0, 1, 2}, {}, {}));
}

TEST(ProbDistributionsUniform, distributionTest) {
  check_quantiles_real_real(UniformTestRig({-1.7, -0.1, 0.0, 0.5},
                                           {-2, -1, 0, 1}, {}, {}, {1.1, 3.8},
                                           {2, 6}, {}, {}));

  check_quantiles_real_real(UniformTestRig({-7.5, -5.1, -4.0}, {-7, -6, -4}, {},
                                           {}, {-2.2, 0.0, 1.0}, {-3, 0, 1, 2},
                                           {}, {}));
}

TEST(ProbDistributionsUniform, error_check) {
  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::uniform_rng(1.0, 2.0, rng));

  EXPECT_THROW(
      stan::math::uniform_rng(stan::math::negative_infinity(), 2.0, rng),
      std::domain_error);
  EXPECT_THROW(stan::math::uniform_rng(1, stan::math::positive_infinity(), rng),
               std::domain_error);
  EXPECT_THROW(stan::math::uniform_rng(2, 1, rng), std::domain_error);
}

TEST(ProbDistributionsUniform, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000;
  int K = stan::math::round(2 * std::pow(N, 0.4));

  // Generate some samples.
  std::vector<double> samples;
  for (int i = 0; i < N; ++i) {
    samples.push_back(stan::math::uniform_rng(1.0, 2.0, rng));
  }

  // Generate quantiles for the uniform distribution.
  std::vector<double> quantiles;
  for (int i = 1; i <= K; ++i) {
    double frac = static_cast<double>(i) / K;
    quantiles.push_back(1.0 + frac);
  }

  // Assert that they match.
  assert_matches_quantiles(samples, quantiles, 1e-6);
}
