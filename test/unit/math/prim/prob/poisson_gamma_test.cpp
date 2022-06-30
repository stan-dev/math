#include <stan/math/prim.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/prob/VectorIntRNGTestRig.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

class PoissonGammaTestRig : public VectorIntRNGTestRig {
 public:
  PoissonGammaTestRig()
      : VectorIntRNGTestRig(10000, 10, {0, 1, 2, 3, 4, 5, 6}, {0.1, 1.7, 3.99},
                            {1, 2, 3}, {-2.1, -0.5, 0.0}, {-3, -1, 0},
                            {0.1, 1.1, 4.99}, {1, 2, 3}, {-3.0, -2.0, 0.0},
                            {-3, -1, 0}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& alpha, const T2& beta, const T3&,
                        T_rng& rng) const {
    return stan::math::poisson_gamma_rng(alpha, beta, rng);
  }

  template <typename T1>
  double pmf(int y, T1 alpha, double beta, double) const {
    return std::exp(stan::math::poisson_gamma_lpmf(y, alpha, beta));
  }
};

TEST(ProbDistributionsPoissonGamma, errorCheck) {
  check_dist_throws_all_types(PoissonGammaTestRig());
}

TEST(ProbDistributionsPoissonGamma, distributionCheck) {
  check_counts_real_real(PoissonGammaTestRig());
}


TEST(ProbDistributionsPoissonGamma, values) {
  using stan::math::poisson_gamma_lpmf;

  // Reference values calculated by extraDistr::dgpois
  EXPECT_FLOAT_EQ(poisson_gamma_lpmf(10, 4, 2), -6.95199150829391);
  EXPECT_FLOAT_EQ(poisson_gamma_lpmf(1, 1, 1), -1.38629436111989);
  EXPECT_FLOAT_EQ(poisson_gamma_lpmf(0, 1, 1), -0.693147180559945);
  EXPECT_FLOAT_EQ(poisson_gamma_lpmf(1000000, 1, 1), -693147.873707126);
}
