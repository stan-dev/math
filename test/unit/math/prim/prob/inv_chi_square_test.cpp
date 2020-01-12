#include <stan/math/prim.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/prob/util.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

class InvChiSquareTestRig : public VectorRealRNGTestRig {
 public:
  InvChiSquareTestRig()
      : VectorRealRNGTestRig(10000, 10, {0.1, 1.0, 2.5, 4.0}, {1, 2, 3, 4},
                             {-2.7, -1.5, -0.5, 0.0}, {-3, -2, -1, 0}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& nu, const T2&, const T3&, T_rng& rng) const {
    return stan::math::inv_chi_square_rng(nu, rng);
  }

  std::vector<double> generate_quantiles(double nu, double, double) const {
    std::vector<double> quantiles;
    double K = stan::math::round(2 * std::pow(N_, 0.4));
    boost::math::inverse_chi_squared_distribution<> dist(nu);

    for (int i = 1; i < K; ++i) {
      double frac = i / K;
      quantiles.push_back(quantile(dist, frac));
    }
    quantiles.push_back(std::numeric_limits<double>::max());

    return quantiles;
  }
};

TEST(ProbDistributionsInvChiSquare, errorCheck) {
  check_dist_throws_all_types(InvChiSquareTestRig());
}

TEST(ProbDistributionsInvChiSquare, distributionTest) {
  check_quantiles_real(InvChiSquareTestRig());
}

TEST(ProbDistributionsInvChiSquare, error_check) {
  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::inv_chi_square_rng(4.0, rng));

  EXPECT_THROW(stan::math::inv_chi_square_rng(-4.0, rng), std::domain_error);
  EXPECT_THROW(
      stan::math::inv_chi_square_rng(stan::math::positive_infinity(), rng),
      std::domain_error);
}

TEST(ProbDistributionsInvChiSquare, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000;
  int K = stan::math::round(2 * std::pow(N, 0.4));

  std::vector<double> samples;
  for (int i = 0; i < N; ++i) {
    samples.push_back(stan::math::inv_chi_square_rng(4.0, rng));
  }

  // Generate quantiles from boost's inverse chi squared distribution
  boost::math::inverse_chi_squared_distribution<> dist(4.0);
  std::vector<double> quantiles;
  for (int i = 1; i < K; ++i) {
    double frac = static_cast<double>(i) / K;
    quantiles.push_back(quantile(dist, frac));
  }
  quantiles.push_back(std::numeric_limits<double>::max());

  // Assert that they match
  assert_matches_quantiles(samples, quantiles, 1e-6);
}
