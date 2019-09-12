#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/prob/vector_rng_test_helper.hpp>
#include <limits>
#include <vector>

class NormalTestRig : public VectorRealRNGTestRig {
 public:
  /*
   * The default NormalTestRig constructor initializes the TestRig with
   * valid and invalid parameters for a random number generator with two
   * arguments.
   */
  NormalTestRig()
      : VectorRealRNGTestRig(
            10000,  // Number of samples used for quantiles tests
            10,     // Length of vectors for vectorization tests
            {-2.5, -1.7, -0.1, 0.0, 2.0, 5.8},  // Valid values for p1
            {-3, -2, -1, 0, 2, 6},              // Valid integer values for p1
            {}, {}, {0.1, 1.0, 2.5, 4.0}, {1, 2, 3, 4}, {-2.7, -1.5, -0.5, 0.0},
            {-3, -2, -1, 0}) {}

  /*
   * This function wraps up the random number generator for testing.
   *
   * The tested rng can have up to three parameters. Any unused parameters can
   * be ignored.
   */
  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& mean, const T2& sd, const T3& unused,
                        T_rng& rng) const {
    return stan::math::normal_rng(mean, sd, rng);
  }

  /*
   * This function builds the quantiles that we will supply to
   * assert_matches_quantiles to test the normal_rng
   */
  std::vector<double> generate_quantiles(double mu, double sigma,
                                         double unused) const {
    std::vector<double> quantiles;
    double K = stan::math::round(2 * std::pow(N_, 0.4));
    boost::math::normal_distribution<> dist(mu, sigma);

    for (int i = 1; i < K; ++i) {
      double frac = i / K;
      quantiles.push_back(quantile(dist, frac));
    }
    quantiles.push_back(std::numeric_limits<double>::max());

    return quantiles;
  }
};

TEST(ProbDistributionsNormal, errorCheck) {
  /*
   * This test verifies that normal_rng throws errors in the right places.
   *
   * It does so by calling test_rig::generate_samples for all possible
   * combinations of calling arguments.
   */
  check_dist_throws_all_types(NormalTestRig());
}

TEST(ProbDistributionsNormal, distributionTest) {
  /*
   * This test checks that the normal_rng is actually generating numbers from
   * the correct distributions. Quantiles are computed from
   * test_rig::generate_quantiles
   *
   * It does so for all possible combinations of calling arguments.
   */
  check_quantiles_real_real(NormalTestRig());
}
