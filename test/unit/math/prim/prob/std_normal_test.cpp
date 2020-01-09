#include <stan/math/prim.hpp>
#include <test/unit/math/prim/prob/util.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

class StdNormalTestRig : public VectorRealRNGTestRig {
 public:
  /*
   * The default StdNormalTestRig constructor initializes the TestRig with
   * valid and invalid parameters for a random number generator with no
   * arguments.
   */
  StdNormalTestRig() : VectorRealRNGTestRig(10000, 10) {}

  /*
   * This function wraps up the random number generator for testing.
   *
   * The tested rng can have up to three parameters. Any unused parameters can
   * be ignored.
   */
  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1&, const T2&, const T3&, T_rng& rng) const {
    return stan::math::std_normal_rng(rng);
  }

  /*
   * This function builds the quantiles that we will supply to
   * assert_matches_quantiles to test the std_normal_rng
   */
  std::vector<double> generate_quantiles(double, double, double) const {
    std::vector<double> quantiles;
    double K = stan::math::round(2 * std::pow(N_, 0.4));
    boost::math::normal_distribution<> dist(0, 1);

    for (int i = 1; i < K; ++i) {
      double frac = i / K;
      quantiles.push_back(quantile(dist, frac));
    }
    quantiles.push_back(std::numeric_limits<double>::max());

    return quantiles;
  }
};

TEST(ProbDistributionsStdNormal, errorCheck) {
  /*
   * This test verifies that std_normal_rng throws errors in the right places.
   *
   * It does so by calling test_rig::generate_samples for all possible
   * combinations of calling arguments.
   */
  check_dist_throws_all_types(StdNormalTestRig());
}

TEST(ProbDistributionsStdNormal, distributionTest) {
  /*
   * This test checks that the std_normal_rng is actually generating numbers
   * from the correct distributions. Quantiles are computed from
   * test_rig::generate_quantiles
   *
   * It does so for all possible combinations of calling arguments.
   */
  check_quantiles_no_params(StdNormalTestRig());
}

TEST(ProbDistributionsStdNormal, error_check) {
  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::std_normal_rng(rng));
}

TEST(ProbDistributionsStdNormal, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000;
  int K = stan::math::round(2 * std::pow(N, 0.4));

  std::vector<double> samples;
  for (int i = 0; i < N; ++i) {
    samples.push_back(stan::math::std_normal_rng(rng));
  }

  // Generate quantiles from boost's normal distribution
  boost::math::normal_distribution<> dist(0.0, 1.0);
  std::vector<double> quantiles;
  for (int i = 1; i < K; ++i) {
    double frac = static_cast<double>(i) / K;
    quantiles.push_back(quantile(dist, frac));
  }
  quantiles.push_back(std::numeric_limits<double>::max());

  // Assert that they match
  assert_matches_quantiles(samples, quantiles, 1e-6);
}
